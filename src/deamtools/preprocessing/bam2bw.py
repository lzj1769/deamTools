from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pyBigWig
import pysam

from deamtools.utils import get_chrom_sizes_from_bam, get_chrom_sizes_from_file

logger = logging.getLogger(__name__)


def _load_regions(bed_path: str) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            regions.setdefault(chrom, []).append((start, end))
    # Merge overlapping intervals to avoid double-counting reads
    for chrom in regions:
        ivs = sorted(regions[chrom])
        merged: list[tuple[int, int]] = [ivs[0]]
        for start, end in ivs[1:]:
            if start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        regions[chrom] = merged
    return regions


def _count_deamination_on_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    chrom_size: int,
    regions: list[tuple[int, int]] | None,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
    mode: str = "count",
) -> tuple[str, np.ndarray]:
    """Per-chromosome deamination counter.

    For ``mode="count"`` (default), returns the number of C->T (forward) /
    G->A (reverse) editing events at each reference base.

    For ``mode="ratio"``, returns the per-base conversion ratio
    ``events / informative_coverage``, where informative coverage is the number
    of reads that contributed a usable C/T base at a reference C (forward
    reads) or a G/A base at a reference G (reverse reads). Bases that are
    neither C nor T (or G nor A on the reverse strand) — sequencing errors,
    SNPs, indels — are excluded from the denominator. Positions with no
    informative coverage produce a ratio of 0.
    """
    if mode not in ("count", "ratio"):
        raise ValueError(f"mode must be 'count' or 'ratio', got {mode!r}")

    events = np.zeros(chrom_size, dtype=np.float32)
    coverage = np.zeros(chrom_size, dtype=np.float32) if mode == "ratio" else None

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        fetch_regions = regions if regions is not None else [(0, chrom_size)]

        for region_start, region_end in fetch_regions:
            region_end = min(region_end, chrom_size)
            ref_seq = fasta.fetch(chrom, region_start, region_end).upper()

            for read in bam.fetch(chrom, region_start, region_end):
                if (
                    read.is_unmapped
                    or read.is_duplicate
                    or read.is_qcfail
                    or read.is_secondary
                    or read.is_supplementary
                ):
                    continue
                if read.mapping_quality < min_mapq:
                    continue

                seq = read.query_sequence
                if seq is None:
                    continue

                quals = read.query_qualities
                is_reverse = read.is_reverse

                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if ref_pos < region_start or ref_pos >= region_end:
                        continue
                    if quals is not None and quals[query_pos] < min_baseq:
                        continue

                    ref_base = ref_seq[ref_pos - region_start]
                    read_base = seq[query_pos]

                    # Forward strand: C->T deamination at reference C.
                    # Reverse strand: G->A (deamination of C on the template strand).
                    if is_reverse:
                        if ref_base == "G":
                            if read_base == "A":
                                events[ref_pos] += 1
                            if coverage is not None and read_base in ("G", "A"):
                                coverage[ref_pos] += 1
                    else:
                        if ref_base == "C":
                            if read_base == "T":
                                events[ref_pos] += 1
                            if coverage is not None and read_base in ("C", "T"):
                                coverage[ref_pos] += 1

    if extend_size > 0:
        kernel = np.ones(2 * extend_size + 1, dtype=np.float32)
        events = np.convolve(events, kernel, mode="same")
        if coverage is not None:
            coverage = np.convolve(coverage, kernel, mode="same")

    if mode == "ratio":
        out = np.zeros_like(events)
        np.divide(events, coverage, out=out, where=coverage > 0)
        return chrom, out

    return chrom, events


def run_bam2bw(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    chrom_sizes_path: str | None = None,
    bed_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    extend_size: int = 0,
    threads: int = 1,
    mode: str = "count",
) -> None:
    if mode not in ("count", "ratio"):
        raise ValueError(f"mode must be 'count' or 'ratio', got {mode!r}")

    logger.info(f"Running bam2bw (mode={mode})")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    if chrom_sizes_path is not None:
        chrom_sizes = get_chrom_sizes_from_file(chrom_sizes_path)
    else:
        logger.info("Inferring chromosome sizes from BAM header")
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            chrom_sizes = get_chrom_sizes_from_bam(bam)

    bed_regions: dict[str, list[tuple[int, int]]] | None = None
    if bed_path is not None:
        logger.info(f"Regions: {bed_path}")
        bed_regions = _load_regions(bed_path)
        logger.info(f"  {sum(len(v) for v in bed_regions.values())} intervals on "
                    f"{len(bed_regions)} chromosome(s)")

    chroms_to_process = [c for c in chrom_sizes if bed_regions is None or c in bed_regions]
    logger.info(f"Processing {len(chroms_to_process)} chromosome(s) "
                f"with {threads} thread(s)")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    results: dict[str, np.ndarray] = {}

    def _process(chrom: str) -> tuple[str, np.ndarray]:
        regions = bed_regions.get(chrom) if bed_regions is not None else None
        return _count_deamination_on_chrom(
            bam_path=bam_path,
            fasta_path=fasta_path,
            chrom=chrom,
            chrom_size=chrom_sizes[chrom],
            regions=regions,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
            mode=mode,
        )

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_process, chrom): chrom for chrom in chroms_to_process}
        for future in as_completed(futures):
            chrom, signal = future.result()
            if mode == "count":
                logger.info(f"  {chrom}: {int(signal.sum())} deamination event(s)")
            else:
                nonzero = int(np.count_nonzero(signal))
                logger.info(f"  {chrom}: {nonzero} position(s) with non-zero ratio")
            results[chrom] = signal

    logger.info(f"Writing {output_path}")
    with pyBigWig.open(output_path, "w") as bw:
        bw.addHeader(list(chrom_sizes.items()))
        for chrom in chroms_to_process:
            signal = results[chrom]
            nonzero = np.nonzero(signal)[0]
            if len(nonzero) == 0:
                continue
            bw.addEntries(
                chrom,
                nonzero.tolist(),
                values=signal[nonzero].tolist(),
                span=1,
            )

    logger.info("Done")
