"""Build a sequence-logo of the deaminase motif from a per-base BigWig.

Consumes a single-base-resolution BigWig of editing-event counts (as
produced by :mod:`deamtools.preprocessing.bam2bw` in count mode) and the
reference FASTA. For each base with a non-zero count the surrounding
window of reference sequence is read from the FASTA and contributed to the
position-weight matrix weighted by the count. The centre base — the
editing site itself — is excluded so the logo reflects the enzyme's
*flanking* sequence preference.

The deamination patterns ``C->T`` (forward strand) and ``G->A`` (reverse
strand) are unified by reverse-complementing the window when the centre
base is ``G``, so the resulting logo is always rendered in the
``C->T`` orientation.

The PWM is converted to information content (bits) and rendered with
:mod:`logomaker`. The bit-score matrix is also written to a sibling
``.csv`` next to the plot.

Inspired by ``plot_motif.py`` from the ACCESS-ATAC-seq project:
https://github.com/pinellolab/ACCESS-ATAC-seq-deprecated/blob/main/plotting/plot_motif.py

Notes
-----
The input BigWig should be produced with ``--extend_size 0`` (the
default). With a non-zero extension each editing event is broadcast to
neighbouring bases, which would inflate the count at non-editing
positions and skew the resulting logo.
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
import pyBigWig  # noqa: E402
import pysam  # noqa: E402

from deamtools.utils import _load_regions

logger = logging.getLogger(__name__)

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    """Return the reverse complement of ``seq`` (A<->T, C<->G, case preserved)."""
    return seq.translate(_RC_TABLE)[::-1]


def _empty_pwm(window_size: int) -> dict[str, list[float]]:
    """Return a zeroed-out per-base count PWM of length ``window_size``."""
    return {b: [0.0] * window_size for b in "ACGTN"}


def _pwm_to_information_df(
    pwm: dict[str, list[float]],
    window_size: int,
) -> pd.DataFrame:
    """Convert a per-base count PWM to bits-per-base, ready for logomaker.

    Information content per position: ``I_i = log2(K) - H_i``, where ``K``
    is the number of base classes (4 after dropping ``N``) and ``H_i`` is
    the Shannon entropy. Each base's bit score at position i is
    ``p_i * I_i``.

    Parameters
    ----------
    pwm : dict[str, list[float]]
        Mapping of base (one of ``"A"``, ``"C"``, ``"G"``, ``"T"``,
        ``"N"``) to per-column counts. ``"N"`` is dropped before
        computing information content.
    window_size : int
        Length of the window. Used to set the centred integer position
        index ``[-window_size // 2, ..., window_size // 2 - 1]`` on the
        returned DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame with ``window_size`` rows and four columns
        ``("A", "C", "G", "T")`` carrying bit-score values; ready for
        :class:`logomaker.Logo`.
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


def _accumulate_pwm_from_bigwig(
    bw: "pyBigWig.pyBigWig",
    fasta: pysam.FastaFile,
    regions_by_chrom: dict[str, list[tuple[int, int]]],
    window_size: int,
) -> dict[str, list[float]]:
    """Accumulate a per-base count PWM from a BigWig of editing-event counts.

    For each region, fetches the per-base values from ``bw`` and the
    reference sequence (with the required ``window_size // 2`` flanking
    margin) from ``fasta``. Every position with a non-zero count is treated
    as an editing site whose surrounding window contributes to the PWM
    weighted by the count. Positions where the reference base is neither
    ``C`` nor ``G`` are skipped (they cannot be deamination sites).

    Parameters
    ----------
    bw : pyBigWig.pyBigWig
        Open BigWig handle (read mode).
    fasta : pysam.FastaFile
        Open reference FASTA handle.
    regions_by_chrom : dict[str, list[tuple[int, int]]]
        Mapping of chromosome name to a list of half-open ``[start, end)``
        intervals to process.
    window_size : int
        Window length in base pairs around each editing site. The centre
        index is ``window_size // 2``; that column is left at zero so the
        editing base itself doesn't dominate the logo.

    Returns
    -------
    dict[str, list[float]]
        Per-base PWM with one entry per ``ACGTN`` key.
    """
    pwm = _empty_pwm(window_size)
    half = window_size // 2

    for chrom, intervals in regions_by_chrom.items():
        try:
            chrom_size = fasta.get_reference_length(chrom)
        except KeyError:
            logger.warning(
                "  %s: not present in FASTA index; skipping", chrom
            )
            continue
        if chrom not in bw.chroms():
            logger.warning(
                "  %s: not present in BigWig header; skipping", chrom
            )
            continue

        for region_start, region_end in intervals:
            region_end = min(region_end, chrom_size)
            if region_end <= region_start:
                continue

            counts = bw.values(chrom, region_start, region_end, numpy=True)

            # Fetch the reference span with enough flanking margin to cover
            # every window centred inside the region.
            fetch_start = max(0, region_start - half)
            fetch_end = min(chrom_size, region_end + half)
            ref_seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()

            nonzero = np.flatnonzero(np.nan_to_num(counts, nan=0.0))
            for offset in nonzero:
                count = float(counts[offset])
                ref_pos = region_start + int(offset)

                centre = ref_seq[ref_pos - fetch_start]
                if centre not in ("C", "G"):
                    continue

                p1 = ref_pos - half
                p2 = p1 + window_size
                if p1 < 0 or p2 > chrom_size:
                    continue

                motif_seq = ref_seq[p1 - fetch_start : p2 - fetch_start]
                if len(motif_seq) != window_size:
                    continue
                if centre == "G":
                    motif_seq = _revcomp(motif_seq)

                for j, b in enumerate(motif_seq):
                    if j == half:
                        continue  # skip the deaminated base itself
                    if b in pwm:
                        pwm[b][j] += count

    return pwm


def _resolve_regions(
    bw: "pyBigWig.pyBigWig",
    bed_path: str | None,
) -> dict[str, list[tuple[int, int]]]:
    """Decide which (chrom, start, end) intervals to process.

    If ``bed_path`` is supplied the BED file is parsed via
    :func:`deamtools.utils._load_regions`. Otherwise one interval per
    chromosome is derived from the BigWig header so the whole genome is
    covered.

    Parameters
    ----------
    bw : pyBigWig.pyBigWig
        Open BigWig handle, used to derive whole-genome bounds when
        ``bed_path`` is not given.
    bed_path : str, optional
        Path to a BED file. When ``None``, the BigWig's chromosomes are
        used.

    Returns
    -------
    dict[str, list[tuple[int, int]]]
        Per-chromosome list of ``[start, end)`` intervals.
    """
    if bed_path is not None:
        df = _load_regions(bed_path)
        regions_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, group in df.groupby("chrom", sort=False):
            regions_by_chrom[str(chrom)] = list(
                zip(group["start"].astype(int), group["end"].astype(int))
            )
        return regions_by_chrom

    return {chrom: [(0, int(size))] for chrom, size in bw.chroms().items()}


def _plot(df: pd.DataFrame, output_path: str) -> None:
    """Render the bit-score DataFrame as a logomaker logo and save to ``output_path``."""
    fig, ax = plt.subplots(1, 1, figsize=(4, 2))
    ax.set_xlabel("Distance from motif center")
    ax.set_ylabel("Bit score")
    logo = logomaker.Logo(df, ax=ax, baseline_width=0)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.xaxis.set_ticks_position("none")
    logo.ax.xaxis.set_tick_params(pad=-1)
    ax.set_title("deamination motif")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_plot_motif(
    bigwig_path: str,
    fasta_path: str,
    output_path: str,
    bed_path: str | None = None,
    window_size: int = 10,
) -> pd.DataFrame:
    """Build and plot the deamination motif logo from a BigWig of edit counts.

    Reads per-base editing-event counts from ``bigwig_path`` (typically
    produced by ``deamtools bam2bw --mode count``), pairs them with the
    reference window from ``fasta_path``, and writes a sequence logo to
    ``output_path``. A sibling ``.csv`` with the bit-score matrix is also
    written.

    Parameters
    ----------
    bigwig_path : str
        Per-base BigWig of editing-event counts.
    fasta_path : str
        Reference FASTA indexed with ``samtools faidx``.
    output_path : str
        Output plot path. Format is inferred from the extension
        (``.png``, ``.pdf``, ``.svg``, ...).
    bed_path : str, optional
        Optional BED file restricting analysis to a subset of intervals.
        When ``None`` the whole genome (per the BigWig header) is used.
    window_size : int, default 10
        Window size in base pairs around each editing site. Must be at
        least 2.

    Returns
    -------
    pandas.DataFrame
        The bit-score matrix that was rendered.

    Raises
    ------
    ValueError
        If ``window_size < 2``.
    RuntimeError
        If no editing sites are observed within the requested intervals
        (so no motif can be built).
    """
    if window_size < 2:
        raise ValueError(f"window_size must be >= 2, got {window_size}")

    logger.info(f"Running plot_motif (window_size={window_size})")
    logger.info(f"BIGWIG: {bigwig_path}")
    logger.info(f"FASTA:  {fasta_path}")

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with pyBigWig.open(bigwig_path, "r") as bw, pysam.FastaFile(fasta_path) as fasta:
        regions_by_chrom = _resolve_regions(bw, bed_path)
        n_intervals = sum(len(v) for v in regions_by_chrom.values())
        logger.info(
            f"Processing {n_intervals} interval(s) on "
            f"{len(regions_by_chrom)} chromosome(s)"
        )

        pwm = _accumulate_pwm_from_bigwig(
            bw=bw,
            fasta=fasta,
            regions_by_chrom=regions_by_chrom,
            window_size=window_size,
        )

    counts_total = sum(pwm["A"]) + sum(pwm["C"]) + sum(pwm["G"]) + sum(pwm["T"])
    if counts_total == 0:
        raise RuntimeError(
            "No editing sites observed; cannot build a motif logo. "
            "Check the BigWig, BED, and that the BigWig is in count mode."
        )

    df = _pwm_to_information_df(pwm, window_size)

    logger.info(f"Writing {output_path}")
    _plot(df, output_path)

    csv_path = os.path.splitext(output_path)[0] + ".csv"
    df.to_csv(csv_path)
    logger.info(f"Writing {csv_path}")

    logger.info("Done")
    return df
