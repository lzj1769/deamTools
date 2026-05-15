"""Build a sequence-logo of the enzyme motif around accessibility events.

Two modes are supported:

  * ``access`` (default) — places the window around each C->T (forward read) or
    G->A (reverse read) deamination site and tallies the surrounding reference
    bases. The deaminated base itself (the window centre) is excluded so the
    logo reflects flanking sequence preference.
  * ``atac`` — places the window around the Tn5 cut site of each read
    (``read_start + 4`` for forward reads, ``read_end - 4`` for reverse reads).

The resulting position-weight matrix is converted to information content
(bits) and rendered via ``logomaker``. A CSV of the bit-score matrix is
written next to the plot.

Inspired by ``plot_motif.py`` from the ACCESS-ATAC-seq project:
https://github.com/pinellolab/ACCESS-ATAC-seq-deprecated/blob/main/plotting/plot_motif.py
"""

from __future__ import annotations

import logging
import os

import matplotlib

matplotlib.use("Agg")  # non-interactive backend; works in headless envs

import logomaker  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pysam  # noqa: E402

from deamtools.utils import _load_regions, get_chrom_sizes_from_bam

logger = logging.getLogger(__name__)

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _empty_pwm(window_size: int) -> dict[str, list[float]]:
    return {b: [0.0] * window_size for b in "ACGTN"}


def _passes_basic_filters(read) -> bool:
    return not (
        read.is_unmapped
        or read.is_duplicate
        or read.is_qcfail
        or read.is_secondary
        or read.is_supplementary
    )


def _pwm_to_information_df(
    pwm: dict[str, list[float]],
    window_size: int,
) -> pd.DataFrame:
    """Convert a per-base count PWM to bits-per-base, ready for logomaker.

    Information content per position: ``I_i = log2(K) - H_i``, where ``K`` is
    the number of base classes (4 after dropping ``N``) and ``H_i`` is the
    Shannon entropy. Each base's bit score at position i is ``p_i * I_i``.
    """
    df = pd.DataFrame(pwm)[["A", "C", "G", "T"]]
    totals = df.sum(axis=1).replace(0, np.nan)
    probs = df.div(totals, axis=0).fillna(0.0)
    p = probs.to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(p > 0, np.log2(p, where=p > 0), 0.0)
    h = -np.sum(p * log_p, axis=1)
    info = np.log2(p.shape[1]) - h
    bits = probs.multiply(info, axis=0)
    bits.index = np.arange(window_size) - window_size // 2
    bits.index.name = "position"
    return bits


def _atac_pwm(
    bam_path: str,
    fasta_path: str,
    regions_by_chrom: dict[str, list[tuple[int, int]]],
    window_size: int,
    min_mapq: int,
) -> dict[str, list[float]]:
    pwm = _empty_pwm(window_size)
    half = window_size // 2

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        chrom_sizes = {c: fasta.get_reference_length(c) for c in regions_by_chrom}
        for chrom, intervals in regions_by_chrom.items():
            chrom_size = chrom_sizes[chrom]
            for region_start, region_end in intervals:
                for read in bam.fetch(chrom, region_start, region_end):
                    if not _passes_basic_filters(read):
                        continue
                    if read.mapping_quality < min_mapq:
                        continue
                    cut_site = (
                        read.reference_end - 4 if read.is_reverse
                        else read.reference_start + 4
                    )
                    p1 = cut_site - half
                    p2 = p1 + window_size
                    if p1 < 0 or p2 > chrom_size:
                        continue
                    seq = fasta.fetch(chrom, p1, p2).upper()
                    if len(seq) != window_size:
                        continue
                    if read.is_reverse:
                        seq = _revcomp(seq)
                    for i, b in enumerate(seq):
                        if b in pwm:
                            pwm[b][i] += 1
    return pwm


def _access_pwm(
    bam_path: str,
    fasta_path: str,
    regions_by_chrom: dict[str, list[tuple[int, int]]],
    window_size: int,
    min_mapq: int,
    min_baseq: int,
) -> dict[str, list[float]]:
    pwm = _empty_pwm(window_size)
    half = window_size // 2

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        chrom_sizes = {c: fasta.get_reference_length(c) for c in regions_by_chrom}
        for chrom, intervals in regions_by_chrom.items():
            chrom_size = chrom_sizes[chrom]
            for region_start, region_end in intervals:
                ref_seq = fasta.fetch(chrom, region_start, region_end).upper()
                for read in bam.fetch(chrom, region_start, region_end):
                    if not _passes_basic_filters(read):
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

                        if is_reverse:
                            if ref_base != "G" or read_base != "A":
                                continue
                        else:
                            if ref_base != "C" or read_base != "T":
                                continue

                        p1 = ref_pos - half
                        p2 = p1 + window_size
                        if p1 < 0 or p2 > chrom_size:
                            continue
                        motif_seq = fasta.fetch(chrom, p1, p2).upper()
                        if len(motif_seq) != window_size:
                            continue
                        if is_reverse:
                            motif_seq = _revcomp(motif_seq)

                        for j, b in enumerate(motif_seq):
                            if j == half:
                                continue  # skip the deaminated base itself
                            if b in pwm:
                                pwm[b][j] += 1
    return pwm


def _resolve_regions(
    bam_path: str,
    bed_path: str | None,
) -> dict[str, list[tuple[int, int]]]:
    if bed_path is not None:
        df = _load_regions(bed_path)
        regions_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, group in df.groupby("chrom", sort=False):
            regions_by_chrom[str(chrom)] = list(
                zip(group["start"].astype(int), group["end"].astype(int))
            )
        return regions_by_chrom

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_sizes = get_chrom_sizes_from_bam(bam)
    return {c: [(0, sz)] for c, sz in chrom_sizes.items()}


def _plot(df: pd.DataFrame, output_path: str, mode: str) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(4, 2))
    ax.set_xlabel("Distance from motif center")
    ax.set_ylabel("Bit score")
    logo = logomaker.Logo(df, ax=ax, baseline_width=0)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.xaxis.set_ticks_position("none")
    logo.ax.xaxis.set_tick_params(pad=-1)
    ax.set_title(f"{mode} motif")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_plot_motif(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    bed_path: str | None = None,
    mode: str = "access",
    window_size: int = 10,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> pd.DataFrame:
    """Build and plot the per-base bit-score logo for accessibility events.

    Returns the bit-score DataFrame that was rendered, after also writing it
    as ``<output_path-without-extension>.csv`` next to the plot.
    """
    if mode not in ("access", "atac"):
        raise ValueError(f"mode must be 'access' or 'atac', got {mode!r}")
    if window_size < 2:
        raise ValueError(f"window_size must be >= 2, got {window_size}")

    logger.info(f"Running plot_motif (mode={mode}, window_size={window_size})")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    regions_by_chrom = _resolve_regions(bam_path, bed_path)
    n_intervals = sum(len(v) for v in regions_by_chrom.values())
    logger.info(
        f"Processing {n_intervals} interval(s) on "
        f"{len(regions_by_chrom)} chromosome(s)"
    )

    if mode == "atac":
        pwm = _atac_pwm(
            bam_path=bam_path,
            fasta_path=fasta_path,
            regions_by_chrom=regions_by_chrom,
            window_size=window_size,
            min_mapq=min_mapq,
        )
    else:
        pwm = _access_pwm(
            bam_path=bam_path,
            fasta_path=fasta_path,
            regions_by_chrom=regions_by_chrom,
            window_size=window_size,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
        )

    counts_total = sum(pwm["A"]) + sum(pwm["C"]) + sum(pwm["G"]) + sum(pwm["T"])
    if counts_total == 0:
        raise RuntimeError(
            "No editing/cut sites observed; cannot build a motif logo. "
            "Check the BAM, BED, and quality thresholds."
        )

    df = _pwm_to_information_df(pwm, window_size)

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    logger.info(f"Writing {output_path}")
    _plot(df, output_path, mode)

    csv_path = os.path.splitext(output_path)[0] + ".csv"
    df.to_csv(csv_path)
    logger.info(f"Writing {csv_path}")

    logger.info("Done")
    return df
