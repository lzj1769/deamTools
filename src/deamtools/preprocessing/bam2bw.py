from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pyBigWig
import pysam

from deamtools.utils import (
    _load_regions,
    get_chrom_sizes_from_bam,
    get_chrom_sizes_from_file,
)

logger = logging.getLogger(__name__)


def _get_edit_count(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    extend_size: int,
    min_mapq: int,
    min_baseq: int,
) -> np.ndarray:
    """Tally per-base deamination edit counts in a genomic region.

    Iterates every primary, non-duplicate aligned read overlapping
    ``[start, end)`` on ``chrom`` and counts deamination events at each
    reference base. An event is recorded at any reference position where
    the aligned read base differs from the reference in either of the two
    deamination patterns: ``C -> T`` or ``G -> A``. The check is
    strand-agnostic: the same mismatch pattern is counted regardless of
    whether the read maps to the forward or reverse strand.

    When ``extend_size > 0`` each event is broadcast symmetrically into a
    window of width ``2 * extend_size + 1`` around the editing site
    (clipped to the region boundaries), matching the behaviour of the
    upstream reference implementation
    (https://github.com/pinellolab/ACCESS-ATAC).

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open BAM handle. Must be coordinate-sorted and indexed.
    fasta : pysam.FastaFile
        Open reference FASTA handle. Must be indexed (``.fai`` present).
    chrom : str
        Chromosome name.
    start, end : int
        Half-open ``[start, end)`` interval on ``chrom``.
    extend_size : int
        Symmetric extension width in base pairs. Set to ``0`` to record one
        unit at the exact editing site.
    min_mapq : int
        Skip reads whose mapping quality is strictly below this value.
    min_baseq : int
        Skip individual read bases whose quality is strictly below this
        value.

    Returns
    -------
    numpy.ndarray
        1-D ``float32`` array of length ``end - start``. Each entry is the
        number of deamination events at that base within the region.

    Notes
    -----
    Reads flagged as unmapped, duplicate, QC-fail, secondary, or
    supplementary are skipped before the MAPQ check. Indels are handled
    naturally by :meth:`pysam.AlignedSegment.get_aligned_pairs` with
    ``matches_only=True``: only matched (M/=/X) bases contribute.
    """
    width = end - start
    signal = np.zeros(width, dtype=np.float32)
    ref_seq = fasta.fetch(chrom, start, end).upper()

    for read in bam.fetch(reference=chrom, start=start, end=end):
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

        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if ref_pos < start or ref_pos >= end:
                continue
            if quals is not None and quals[query_pos] < min_baseq:
                continue

            ref_base = ref_seq[ref_pos - start]
            read_base = seq[query_pos]

            if (ref_base == "C" and read_base == "T") or (
                ref_base == "G" and read_base == "A"
            ):
                if extend_size > 0:
                    lo = max(0, ref_pos - start - extend_size)
                    hi = min(ref_pos - start + extend_size + 1, width)
                    signal[lo:hi] += 1
                else:
                    signal[ref_pos - start] += 1

    return signal


def _get_total_coverage(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    min_baseq: int,
) -> np.ndarray:
    """Per-base total read coverage in a genomic region.

    Uses :meth:`pysam.AlignmentFile.count_coverage` to sum the per-base
    counts of A, C, G, and T across all reads overlapping ``[start, end)``.
    Reads flagged as unmapped, duplicate, QC-fail, secondary, or
    supplementary are excluded by pysam's default filter mask.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open BAM handle.
    chrom : str
        Chromosome name.
    start, end : int
        Half-open ``[start, end)`` interval on ``chrom``.
    min_baseq : int
        Per-base quality threshold forwarded to ``count_coverage`` (bases
        with quality strictly below this value are not counted).

    Returns
    -------
    numpy.ndarray
        1-D ``float32`` array of length ``end - start`` giving the total
        ACGT read coverage at each base in the region.
    """
    width = end - start
    a, c, g, t = bam.count_coverage(
        chrom, start, end, quality_threshold=min_baseq
    )
    cov = np.array(a, dtype=np.float32)
    cov += np.array(c, dtype=np.float32)
    cov += np.array(g, dtype=np.float32)
    cov += np.array(t, dtype=np.float32)
    assert cov.shape[0] == width
    return cov


def _signal_for_region(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    start: int,
    end: int,
    mode: str,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
    min_coverage: int,
) -> tuple[str, int, int, np.ndarray]:
    """Compute the per-base signal for a single genomic region.

    Wraps :func:`_get_edit_count` (and, in ratio mode,
    :func:`_get_total_coverage`) so the result can be dispatched to a
    worker thread. The returned tuple includes the input coordinates so
    the orchestrator can assemble outputs in BigWig-sorted order
    independent of completion order.

    Parameters
    ----------
    bam_path, fasta_path : str
        Paths to the BAM and reference FASTA. Both are opened locally so
        the function is safe to call from worker threads (pysam handles
        are not thread-safe).
    chrom : str
        Chromosome name.
    start, end : int
        Half-open ``[start, end)`` interval on ``chrom``.
    mode : {"count", "ratio"}
        Signal type. ``"count"`` returns the edit count (possibly extended
        when ``extend_size > 0``). ``"ratio"`` returns
        ``edit_count / total_coverage`` with positions whose coverage is
        strictly below ``min_coverage`` masked to ``0``.
    min_mapq, min_baseq, extend_size, min_coverage : int
        See :func:`run_bam2bw`.

    Returns
    -------
    chrom : str
        Echoed input chromosome.
    start, end : int
        Echoed input coordinates.
    signal : numpy.ndarray
        1-D ``float32`` array of length ``end - start``.
    """
    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        edits = _get_edit_count(
            bam=bam,
            fasta=fasta,
            chrom=chrom,
            start=start,
            end=end,
            extend_size=extend_size if mode == "count" else 0,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
        )

        if mode == "count":
            return chrom, start, end, edits

        coverage = _get_total_coverage(
            bam=bam,
            chrom=chrom,
            start=start,
            end=end,
            min_baseq=min_baseq,
        )
        coverage = np.where(coverage < min_coverage, 0.0, coverage)

        signal = np.zeros_like(edits)
        np.divide(edits, coverage, out=signal, where=coverage > 0)
        return chrom, start, end, signal


def run_bam2bw(
    bam_path: str,
    fasta_path: str,
    out_dir: str,
    out_name: str,
    chrom_sizes_path: str | None = None,
    bed_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    extend_size: int = 0,
    threads: int = 1,
    mode: str = "count",
    min_coverage: int = 1,
) -> None:
    """Convert a BAM file to a per-base BigWig track of deamination signal.

    Drives the end-to-end pipeline: enumerates genomic regions (either the
    intervals in ``bed_path`` or one region per chromosome derived from the
    BAM header), computes the per-base signal for each region in parallel
    via :func:`_signal_for_region`, and writes a BigWig with one entry per
    non-zero base in sorted order.

    The algorithm mirrors the upstream ACCESS-ATAC reference
    (https://github.com/pinellolab/ACCESS-ATAC): edits are counted in a
    strand-agnostic fashion (``C->T`` or ``G->A`` reference mismatches);
    ``--extend_size`` applies only in count mode; the fraction mode
    denominator is the total ACGT coverage with positions below
    ``min_coverage`` masked to zero.

    Parameters
    ----------
    bam_path : str
        Path to the coordinate-sorted, indexed BAM file.
    fasta_path : str
        Path to the indexed reference FASTA file
        (``samtools faidx``-style ``.fai`` required).
    out_dir : str
        Output directory. Created if it does not already exist.
    out_name : str
        Base name (without extension) for the output; the BigWig is written to
        ``<out_dir>/<out_name>.bw``.
    chrom_sizes_path : str, optional
        Tab-delimited chromosome-sizes file in UCSC format
        (``chrom`` ``<TAB>`` ``size`` per line). When ``None`` (default),
        chromosome sizes are inferred from the BAM header.
    bed_path : str, optional
        BED file restricting analysis to a subset of intervals. When ``None``
        (default), the entire genome is processed. The BED is parsed with
        :func:`_load_regions`, which merges overlapping/adjacent intervals
        to prevent double-counting.
    min_mapq : int, default 20
        Minimum read mapping quality.
    min_baseq : int, default 20
        Minimum base quality at a position to count an editing event.
    extend_size : int, default 0
        Symmetric extension width in base pairs. Only applied in count
        mode: each editing event is broadcast into a window of width
        ``2 * extend_size + 1`` centred on the editing site (clipped to
        the enclosing region). Ignored in ratio mode.
    threads : int, default 1
        Number of worker threads used to process regions in parallel.
    mode : {"count", "ratio"}, default "count"
        Signal to write to the BigWig.

        * ``"count"`` — raw per-base deamination edit count.
        * ``"ratio"`` — per-base conversion ratio
          ``edit_count / total_coverage``. Positions whose total ACGT
          coverage is strictly below ``min_coverage`` are written as
          ``0``.
    min_coverage : int, default 1
        Coverage threshold for ratio mode (ignored when ``mode="count"``).
        Positions whose total ACGT coverage is strictly below this value
        report a ratio of ``0`` rather than a noisy small-denominator
        fraction.

    Returns
    -------
    None
        The BigWig is written to ``<out_dir>/<out_name>.bw`` as a side effect. The
        header lists every chromosome from ``chrom_sizes_path`` (or the BAM
        header); chromosomes with no signal are present in the header but
        carry no entries.

    Raises
    ------
    ValueError
        If ``mode`` is not ``"count"`` or ``"ratio"``.
    FileNotFoundError
        Re-raised from :mod:`pysam` / :mod:`pyBigWig` if any required input
        or index file is missing.

    See Also
    --------
    _count_deamination_on_chrom : Per-chromosome counter used internally.
    _load_regions : BED loader used to restrict processing.
    """
    if mode not in ("count", "ratio"):
        raise ValueError(f"mode must be 'count' or 'ratio', got {mode!r}")

    output_path = os.path.join(out_dir, f"{out_name}.bw")

    logger.info(f"Running bam2bw (mode={mode})")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    if chrom_sizes_path is not None:
        chrom_sizes = get_chrom_sizes_from_file(chrom_sizes_path)
    else:
        logger.info("Inferring chromosome sizes from BAM header")
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            chrom_sizes = get_chrom_sizes_from_bam(bam)

    # Build the (chrom, start, end) work list.
    if bed_path is not None:
        logger.info(f"Regions: {bed_path}")
        bed_df = _load_regions(bed_path)
        bed_df = bed_df[bed_df["chrom"].isin(chrom_sizes)].reset_index(drop=True)
        regions: list[tuple[str, int, int]] = [
            (str(c), int(s), min(int(e), chrom_sizes[c]))
            for c, s, e in zip(bed_df["chrom"], bed_df["start"], bed_df["end"])
        ]
        logger.info(
            f"  {len(regions)} interval(s) on "
            f"{bed_df['chrom'].nunique()} chromosome(s)"
        )
    else:
        logger.info("No BED supplied; using whole genome")
        regions = [(c, 0, chrom_sizes[c]) for c in chrom_sizes]

    # Sort to satisfy pyBigWig's requirement that entries are added in
    # (chrom, start, end) order. We follow the chromosome order in the BAM
    # header, which is what pyBigWig will use for the header itself.
    chrom_order = {c: i for i, c in enumerate(chrom_sizes)}
    regions.sort(key=lambda r: (chrom_order[r[0]], r[1], r[2]))

    logger.info(
        f"Processing {len(regions)} region(s) with {threads} thread(s)"
    )

    os.makedirs(out_dir, exist_ok=True)

    results: dict[tuple[str, int, int], np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(
                _signal_for_region,
                bam_path=bam_path,
                fasta_path=fasta_path,
                chrom=c,
                start=s,
                end=e,
                mode=mode,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                extend_size=extend_size,
                min_coverage=min_coverage,
            )
            for c, s, e in regions
        ]
        for future in as_completed(futures):
            chrom, start, end, signal = future.result()
            results[(chrom, start, end)] = signal

    if mode == "count":
        total = int(sum(int(s.sum()) for s in results.values()))
        logger.info(f"  total deamination event(s): {total}")
    else:
        nonzero = int(sum(int(np.count_nonzero(s)) for s in results.values()))
        logger.info(f"  total position(s) with non-zero ratio: {nonzero}")

    logger.info(f"Writing {output_path}")
    with pyBigWig.open(output_path, "w") as bw:
        bw.addHeader(list(chrom_sizes.items()))
        for chrom, start, end in regions:
            signal = results[(chrom, start, end)]
            nonzero_idx = np.nonzero(signal)[0]
            if len(nonzero_idx) == 0:
                continue
            bw.addEntries(
                chrom,
                (nonzero_idx + start).tolist(),
                values=signal[nonzero_idx].astype(float).tolist(),
                span=1,
            )

    logger.info("Done")
