"""Transcription-factor footprint scoring from a per-base editing BigWig.

For each motif-predicted binding site (a region in the input BED) of width
``L = end - start``, the per-base signal is read over the ``3 * L`` window
``[start - L, end + L)`` and split into three equal parts — the left flank, the
motif centre, and the right flank — and scored as

    fp_score = mean(left flank) + mean(right flank) - mean(centre)

A bound TF shields its motif from editing, so a footprint shows depletion in
the centre relative to the flanks and yields a positive score. Significance is
assessed by permuting the per-base signal within the window many times to build
a null distribution of footprint scores; the p-value is the fraction of
permutations scoring at least as high as the observed score.

The algorithm follows ``05_compute_fp_score.ipynb``. The output is a 6-column,
BED-like table: ``chrom  start  end  name  fp_score  p_value``.
"""

from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pyBigWig

logger = logging.getLogger(__name__)


def _read_bed(regions_path: str) -> list[tuple[str, int, int, str]]:
    """Read a BED of motif sites as ``(chrom, start, end, name)`` (no merging)."""
    records: list[tuple[str, int, int, str]] = []
    with open(regions_path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else "."
            records.append((chrom, start, end, name))
    return records


def _footprint_score(signal: np.ndarray, motif_length: int) -> float:
    """flank means minus centre mean for a ``3 * motif_length`` signal window."""
    left = signal[:motif_length].mean()
    centre = signal[motif_length : 2 * motif_length].mean()
    right = signal[2 * motif_length :].mean()
    return float(left + right - centre)


def _score_chrom(
    bigwig_path: str,
    records: list[tuple[str, int, int, str]],
    n_shuffles: int,
    seed: int | None,
) -> tuple[list[str], int]:
    """Score every record on one chromosome; returns BED rows and a skip count.

    Opens its own BigWig handle so it is safe to call from a worker thread.
    """
    rng = np.random.default_rng(seed)
    rows: list[str] = []
    skipped = 0
    with pyBigWig.open(bigwig_path) as bw:
        chroms = bw.chroms()
        for chrom, start, end, name in records:
            motif_length = end - start
            lo = start - motif_length
            hi = end + motif_length
            chrom_len = chroms.get(chrom)
            # Need the full 3x-motif window inside the chromosome.
            if motif_length <= 0 or chrom_len is None or lo < 0 or hi > chrom_len:
                skipped += 1
                continue

            signal = np.nan_to_num(
                np.asarray(bw.values(chrom, lo, hi), dtype=float)
            )
            fp_score = _footprint_score(signal, motif_length)

            if fp_score <= 0:
                p_value = 1.0
            else:
                # Vectorised null: independently permute the window n_shuffles
                # times and recompute the footprint score for each.
                tiled = np.tile(signal, (n_shuffles, 1))
                shuffled = rng.permuted(tiled, axis=1)
                left = shuffled[:, :motif_length].mean(axis=1)
                centre = shuffled[:, motif_length : 2 * motif_length].mean(axis=1)
                right = shuffled[:, 2 * motif_length :].mean(axis=1)
                bg = left + right - centre
                p_value = float((np.sum(bg >= fp_score) + 1) / (n_shuffles + 1))

            rows.append(
                f"{chrom}\t{start}\t{end}\t{name}\t{fp_score:.6g}\t{p_value:.6g}"
            )
    return rows, skipped


def run_footprint(
    bigwig_path: str,
    regions_path: str,
    out_dir: str,
    out_name: str,
    n_shuffles: int = 1000,
    threads: int = 1,
    seed: int | None = None,
) -> None:
    """Compute footprint scores for motif sites and write a BED-like table.

    Parameters
    ----------
    bigwig_path : str
        Per-base editing BigWig, e.g. produced by ``deamtools bam2bw``.
    regions_path : str
        BED of motif-predicted binding sites (e.g. from ``deamtools match``);
        column 4, if present, is carried through as the site name.
    out_dir : str
        Output directory. Created if it does not exist.
    out_name : str
        Base name (without extension) for the output; writes
        ``<out_dir>/<out_name>.bed`` with columns
        ``chrom  start  end  name  fp_score  p_value``.
    n_shuffles : int, default 1000
        Number of within-window permutations used to build the null
        distribution (only computed for sites with a positive score).
    threads : int, default 1
        Number of worker threads; chromosomes are scored in parallel.
    seed : int, optional
        Base RNG seed for reproducible p-values (each chromosome gets a derived
        seed). When ``None``, results are non-deterministic.
    """
    logger.info("Running footprint")
    logger.info(f"BigWig:  {bigwig_path}")
    logger.info(f"Regions: {regions_path}")

    records = _read_bed(regions_path)
    logger.info(f"  {len(records)} region(s); {n_shuffles} shuffle(s)")

    # Group by chromosome, preserving first-seen order for stable output.
    by_chrom: dict[str, list[tuple[str, int, int, str]]] = {}
    for rec in records:
        by_chrom.setdefault(rec[0], []).append(rec)

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, f"{out_name}.bed")

    results: dict[str, list[str]] = {}
    total_skipped = 0
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {}
        for i, (chrom, recs) in enumerate(by_chrom.items()):
            # Derive a per-chromosome seed so output is reproducible yet the
            # chromosomes use independent random streams.
            chrom_seed = None if seed is None else seed + i
            futures[
                pool.submit(_score_chrom, bigwig_path, recs, n_shuffles, chrom_seed)
            ] = chrom
        for future in as_completed(futures):
            rows, skipped = future.result()
            results[futures[future]] = rows
            total_skipped += skipped

    logger.info(f"Writing {output_path}")
    n_written = 0
    with open(output_path, "w") as out:
        for chrom in by_chrom:
            for row in results.get(chrom, []):
                out.write(row + "\n")
                n_written += 1

    if total_skipped:
        logger.info(
            f"  skipped {total_skipped} region(s) out of bounds / zero-length"
        )
    logger.info(f"  wrote {n_written} footprint(s)")
    logger.info("Done")
