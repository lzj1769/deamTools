from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Sequence
from functools import partial

import numpy as np
import pyBigWig
import pysam

from deamtools.utils import (
    _load_regions,
    get_chrom_sizes_from_bam,
    get_chrom_sizes_from_file,
    iter_fragments,
    merge_fragment_bases,
    run_jobs,
)

logger = logging.getLogger(__name__)


_ACGT = frozenset("ACGT")

# What the track counts: deamination edits, or Tn5 insertion (cut) sites.
EVENT_EDIT = "edit"
EVENT_TN5 = "tn5"
EVENTS = (EVENT_EDIT, EVENT_TN5)

# Tn5 inserts as a dimer that nicks the two strands 9 bp apart, so both
# fragments flanking one insertion carry the same 9-bp duplication [p, p + 9).
# A forward read of the right-hand fragment starts at p; a reverse read of the
# left-hand fragment ends (exclusive) at p + 9. Moving each 4 bp inward from its
# 5' end -- +4 on reference_start, -5 on the exclusive reference_end -- puts
# both on p + 4, the centre base of the duplication: the usual +4/-5 shift.
TN5_SHIFT_FORWARD = 4
TN5_SHIFT_REVERSE = 5


def _passes_filters(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    """Primary, non-duplicate, mapping-quality-passing read?"""
    if (
        read.is_unmapped
        or read.is_duplicate
        or read.is_qcfail
        or read.is_secondary
        or read.is_supplementary
    ):
        return False
    return read.mapping_quality >= min_mapq


def _tn5_cut_site(read: pysam.AlignedSegment) -> int | None:
    """0-based reference position of the Tn5 insertion that made this read.

    That is the read's **5' end**, shifted to the centre of the 9-bp target-site
    duplication: ``reference_start + 4`` for a forward read, and
    ``reference_end - 5`` for a reverse one, whose 5' end is its right-most
    aligned base. See ``TN5_SHIFT_FORWARD`` for the geometry.
    """
    if read.is_reverse:
        end = read.reference_end
        return None if end is None else end - TN5_SHIFT_REVERSE
    return read.reference_start + TN5_SHIFT_FORWARD


def _get_cut_count(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    extend_size: int,
    min_mapq: int,
) -> np.ndarray:
    """Tally per-base Tn5 cut sites in ``[start, end)``.

    Every passing read contributes the cut at its 5' end (:func:`_tn5_cut_site`).
    For paired-end data that is exactly both ends of each fragment -- the two
    mates' 5' ends are the fragment's two ends, one Tn5 insertion each -- and for
    single-end data it is the one end the read was sequenced from. So one rule
    gives both layouts, including a mate whose partner was filtered out, which
    keeps its own real cut.

    Unlike edits, a cut needs no reference base and no base quality: it is a
    property of where the read's alignment starts. Filters and ``extend_size``
    work as they do for edits.
    """
    width = end - start
    signal = np.zeros(width, dtype=np.float32)
    for read in bam.fetch(reference=chrom, start=start, end=end):
        if not _passes_filters(read, min_mapq):
            continue
        cut = _tn5_cut_site(read)
        if cut is None or not start <= cut < end:
            continue
        idx = cut - start
        if extend_size > 0:
            signal[max(0, idx - extend_size) : min(idx + extend_size + 1, width)] += 1
        else:
            signal[idx] += 1
    return signal


def _get_edit_count(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    extend_size: int,
    min_mapq: int,
    min_baseq: int,
    want_coverage: bool = False,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Tally per-base deamination edit counts, and coverage, in a region.

    Iterates every primary, non-duplicate aligned read overlapping
    ``[start, end)`` on ``chrom`` and counts deamination events at each
    reference base. An event is recorded at any reference position where
    the aligned read base differs from the reference in either of the two
    deamination patterns: ``C -> T`` or ``G -> A``. The check is
    strand-agnostic: the same mismatch pattern is counted regardless of
    whether the read maps to the forward or reverse strand.

    Counting is per **fragment**, not per record. Records are grouped with
    :func:`~deamtools.utils.iter_fragments` and each fragment's mates are
    collapsed by :func:`~deamtools.utils.merge_fragment_bases`, so a
    reference position that both mates cover contributes one event and one
    unit of coverage rather than two -- mates overlap whenever the insert is
    shorter than twice the read length, and that overlap is the middle of
    the fragment, not a random subset of positions.

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
    want_coverage : bool, default False
        Also return per-base fragment coverage (positions where the merged
        fragment called an A, C, G or T). Computed in the same pass, so the
        numerator and denominator of a ratio always see the same fragments;
        skipped when not needed, since it is a write per aligned base.

    Returns
    -------
    signal : numpy.ndarray
        1-D ``float32`` array of length ``end - start``. Each entry is the
        number of deamination events at that base within the region.
    coverage : numpy.ndarray or None
        Same shape, giving per-base ACGT fragment coverage; ``None`` unless
        ``want_coverage``.

    Notes
    -----
    Reads flagged as unmapped, duplicate, QC-fail, secondary, or
    supplementary are skipped before the MAPQ check. Indels are handled
    naturally by :meth:`pysam.AlignedSegment.get_aligned_pairs` with
    ``matches_only=True``: only matched (M/=/X) bases contribute.
    """
    width = end - start
    signal = np.zeros(width, dtype=np.float32)
    coverage = np.zeros(width, dtype=np.float32) if want_coverage else None
    ref_seq = fasta.fetch(chrom, start, end).upper()

    def passing_reads() -> Iterator[pysam.AlignedSegment]:
        for read in bam.fetch(reference=chrom, start=start, end=end):
            if _passes_filters(read, min_mapq):
                yield read

    for fragment in iter_fragments(passing_reads()):
        bases = merge_fragment_bases(fragment, min_baseq, start, end)
        for ref_pos, read_base in bases.items():
            idx = ref_pos - start

            if coverage is not None and read_base in _ACGT:
                coverage[idx] += 1

            ref_base = ref_seq[idx]
            if (ref_base == "C" and read_base == "T") or (
                ref_base == "G" and read_base == "A"
            ):
                if extend_size > 0:
                    lo = max(0, idx - extend_size)
                    hi = min(idx + extend_size + 1, width)
                    signal[lo:hi] += 1
                else:
                    signal[idx] += 1

    return signal, coverage


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
    event: str = EVENT_EDIT,
) -> tuple[str, int, int, np.ndarray]:
    """Compute the per-base signal for a single genomic region.

    Wraps :func:`_get_edit_count` so the result can be dispatched to a
    worker process. The returned tuple includes the input coordinates so
    the orchestrator can assemble outputs in BigWig-sorted order
    independent of completion order.

    Parameters
    ----------
    bam_path, fasta_path : str
        Paths to the BAM and reference FASTA. Both are opened locally so
        the function is safe to run in a worker process (pysam handles
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
    [(_, _, _, signal)] = _signal_for_batch(
        bam_path=bam_path,
        fasta_path=fasta_path,
        regions=[(chrom, start, end)],
        mode=mode,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        extend_size=extend_size,
        min_coverage=min_coverage,
        event=event,
    )
    return chrom, start, end, signal


def _region_signal(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile | None,
    chrom: str,
    start: int,
    end: int,
    mode: str,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
    min_coverage: int,
    event: str,
) -> np.ndarray:
    """One region's per-base signal, from handles the caller already opened."""
    if event == EVENT_TN5:
        return _get_cut_count(bam, chrom, start, end, extend_size, min_mapq)

    assert fasta is not None  # opened for every edit-mode batch
    edits, coverage = _get_edit_count(
        bam=bam,
        fasta=fasta,
        chrom=chrom,
        start=start,
        end=end,
        extend_size=extend_size if mode == "count" else 0,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        want_coverage=mode == "ratio",
    )
    if mode == "count":
        return edits

    assert coverage is not None  # set whenever want_coverage was requested
    coverage = np.where(coverage < min_coverage, 0.0, coverage)
    signal = np.zeros_like(edits)
    np.divide(edits, coverage, out=signal, where=coverage > 0)
    return signal


def _signal_for_batch(
    bam_path: str,
    fasta_path: str,
    regions: Sequence[tuple[str, int, int]],
    mode: str,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
    min_coverage: int,
    event: str = EVENT_EDIT,
) -> list[tuple[str, int, int, np.ndarray]]:
    """Compute a batch of regions with one BAM (and FASTA) handle between them.

    This is the unit of work handed to a worker. Opening a BAM reads its whole
    index, so opening one per region -- as bam2bw used to -- made a peak run
    mostly index loading: over 57k HepG2 peaks, counting Tn5 cuts on the
    single-end BAM took as long as counting edits (43 s vs 44 s), though a cut
    is trivial to compute next to an edit. Batches are runs of neighbouring
    regions (see :func:`_batch_regions`), so successive fetches also land in
    nearby parts of the file. Tn5 batches never open the FASTA.
    """
    fasta = None if event == EVENT_TN5 else pysam.FastaFile(fasta_path)
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            return [
                (
                    chrom,
                    start,
                    end,
                    _region_signal(
                        bam,
                        fasta,
                        chrom,
                        start,
                        end,
                        mode,
                        min_mapq,
                        min_baseq,
                        extend_size,
                        min_coverage,
                        event,
                    ),
                )
                for chrom, start, end in regions
            ]
    finally:
        if fasta is not None:
            fasta.close()


# A batch is capped at this many reference bases. On a whole-genome run every
# region is a chromosome, so the cap keeps them one per batch and a worker never
# holds several chromosome-length arrays at once; small peaks pack far below it.
_MAX_BATCH_SPAN = 50_000_000
# Batches per worker, so a slow batch does not leave the other workers idle.
_BATCHES_PER_WORKER = 8


def _batch_regions(
    regions: Sequence[tuple[str, int, int]], workers: int
) -> list[list[tuple[str, int, int]]]:
    """Split sorted regions into contiguous runs of roughly equal total span.

    With one worker there is nothing to balance, so everything goes in as few
    batches as the span cap allows -- the serial path gains as much from
    opening the BAM once as the parallel one does. Order is preserved, so
    concatenating the batches gives ``regions`` back.
    """
    if not regions:
        return []
    total = sum(end - start for _, start, end in regions)
    n_target = 1 if workers <= 1 else workers * _BATCHES_PER_WORKER
    span_goal = max(1, min(_MAX_BATCH_SPAN, -(-total // n_target)))
    batches: list[list[tuple[str, int, int]]] = []
    current: list[tuple[str, int, int]] = []
    span = 0
    for region in regions:
        width = region[2] - region[1]
        if current and span + width > span_goal:
            batches.append(current)
            current, span = [], 0
        current.append(region)
        span += width
    batches.append(current)
    return batches


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
    normalize: bool = False,
    scale_factor: float = 1_000_000.0,
    event: str = EVENT_EDIT,
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

    Both the edit count and the ratio denominator are per **fragment**: the
    mates of a pair are merged before counting, so a reference position that
    both mates cover contributes once, not twice. Numerator and denominator
    are computed in the same pass over the same fragments, and both honour
    ``min_mapq`` -- before 2026-09-10 the denominator came from
    ``count_coverage``, which applies no mapping-quality filter, so ratio mode
    was dividing edits from MAPQ-passing reads by coverage that included reads
    the numerator had excluded.

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
        Number of worker processes used to process regions in parallel.
        With 1 everything runs in this process.
    mode : {"count", "ratio"}, default "count"
        Signal to write to the BigWig.

        * ``"count"`` — raw per-base deamination edit count.
        * ``"ratio"`` — per-base conversion ratio
          ``edit_count / total_coverage``. Positions whose total ACGT
          coverage is strictly below ``min_coverage`` are written as
          ``0``. Edits only; rejected with ``event="tn5"``.
    event : {"edit", "tn5"}, default "edit"
        What the track counts.

        * ``"edit"`` — deamination events (C->T or G->A reference
          mismatches), per fragment, as described above.
        * ``"tn5"`` — Tn5 insertion sites: the 5' end of every passing read,
          shifted +4 (forward) or -5 (reverse, from the exclusive end) onto
          the centre of the 9-bp target-site duplication. A paired-end
          fragment therefore contributes both of its ends and a single-end
          read its start. ``min_baseq`` does not apply; ``min_mapq``,
          the flag filters, ``extend_size`` and ``normalize`` do.
    min_coverage : int, default 1
        Coverage threshold for ratio mode (ignored when ``mode="count"``).
        Positions whose total ACGT coverage is strictly below this value
        report a ratio of ``0`` rather than a noisy small-denominator
        fraction.
    normalize : bool, default False
        Apply reads-per-million-style normalization to the **count**-mode
        signal: every value is scaled by ``scale_factor / total``, where
        ``total`` is the genome-wide sum of the count signal. The written
        track therefore sums to ``scale_factor`` (with ``extend_size=0`` this
        is counts-per-``scale_factor`` edits). Ignored in ratio mode.
    scale_factor : float, default 1_000_000
        Target total for ``--normalize`` (1e6 gives reads/counts-per-million).

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
    _get_edit_count : Per-region fragment counter used internally.
    _load_regions : BED loader used to restrict processing.
    """
    if mode not in ("count", "ratio"):
        raise ValueError(f"mode must be 'count' or 'ratio', got {mode!r}")
    if event not in EVENTS:
        raise ValueError(f"event must be one of {', '.join(EVENTS)}, got {event!r}")
    if event == EVENT_TN5 and mode == "ratio":
        raise ValueError(
            "mode='ratio' divides edits by coverage and has no meaning for Tn5 "
            "cut sites; use mode='count' with event='tn5'"
        )

    output_path = os.path.join(out_dir, f"{out_name}.bw")

    logger.info(f"Running bam2bw (event={event}, mode={mode})")
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
            for c, s, e in zip(
                bed_df["chrom"], bed_df["start"], bed_df["end"], strict=True
            )
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

    batches = _batch_regions(regions, threads)
    logger.info(
        f"Processing {len(regions)} region(s) in {len(batches)} batch(es) "
        f"with {threads} worker(s)"
    )

    os.makedirs(out_dir, exist_ok=True)

    results: dict[tuple[str, int, int], np.ndarray] = {}
    jobs = [
        partial(
            _signal_for_batch,
            bam_path=bam_path,
            fasta_path=fasta_path,
            regions=batch,
            mode=mode,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
            min_coverage=min_coverage,
            event=event,
        )
        for batch in batches
    ]
    for batch_result in run_jobs(jobs, threads):
        for chrom, start, end, signal in batch_result:
            results[(chrom, start, end)] = signal

    norm_factor = 1.0
    if mode == "count":
        total = int(sum(int(s.sum()) for s in results.values()))
        what = "Tn5 cut site(s)" if event == EVENT_TN5 else "deamination event(s)"
        logger.info(f"  total {what}: {total}")
        if normalize:
            norm_factor = scale_factor / total if total > 0 else 0.0
            logger.info(
                f"  normalizing by {total} -> scale_factor {scale_factor:g} "
                f"(factor {norm_factor:g})"
            )
    else:
        if normalize:
            logger.warning("  --normalize is ignored in ratio mode")
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
            values = signal[nonzero_idx].astype(float)
            if norm_factor != 1.0:
                values = values * norm_factor
            bw.addEntries(
                chrom,
                (nonzero_idx + start).tolist(),
                values=values.tolist(),
                span=1,
            )

    logger.info("Done")
