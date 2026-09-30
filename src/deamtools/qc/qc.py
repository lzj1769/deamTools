"""Quality-control metrics for deaminase-based chromatin accessibility data.

Summarises a coordinate-sorted BAM together with its reference FASTA into the
metrics most useful for judging a deaminase footprinting experiment:

* **Read statistics** — totals plus the fraction of duplicate, properly-paired,
  secondary and supplementary reads.
* **Editing statistics** — the genome-wide deamination rate (edits divided by
  the number of editable C/G *opportunities* covered by passing fragments) and
  the distribution of edits per fragment. A high, accessibility-driven edit rate
  is the primary signal that the deaminase treatment worked.

  Editing is counted **per fragment, not per record**: the two mates of a pair
  are merged first (:func:`~deamtools.utils.merge_fragment_bases`), so a reference position both mates
  cover is one observation rather than two. Where the mates disagree — which at
  an overlap means a sequencing error, since library prep turns a deaminated C
  into a real T:A pair that both mates then carry — the higher base quality
  wins. Single-end reads are fragments of one record, so nothing changes there.
* **Trinucleotide context bias** — the edit fraction broken down by the
  trinucleotide centred on the edited cytosine. Edits are called
  strand-agnostically (matching :mod:`deamtools.preprocessing.bam2bw`): any
  reference ``C->T`` or ``G->A`` mismatch counts, and ``G``-centred contexts
  are reverse-complemented so both are reported in the unified ``C``-centred
  orientation. This is
  the enzyme's sequence-preference fingerprint (e.g. DddA's ``TC`` preference).
* **Fragment-length distribution** — from the template length of properly-paired
  read pairs. Paired-end only: the library layout is detected from the head of
  the BAM (:func:`_detect_layout`), and for a single-end library this block and
  the proper-pair metrics are omitted rather than reported as zero.
* **Deaminase sequence motif** — a sequence logo of the reference bases flanking
  edited cytosines, built directly from the editing events (centre excluded,
  ``G->A`` events reverse-complemented to the ``C->T`` orientation). Shows the
  enzyme's flanking-sequence preference.
* **TSS enrichment** *(optional)* — enrichment of insertion sites around
  transcription start sites, computed when a TSS BED is supplied, following the
  ENCODE ATAC-seq pipeline definition (see :func:`_tss_enrichment`). The
  aggregate profile behind the plot is also written as CSV.

Results are written as machine-readable JSON plus a self-contained,
MultiQC-style HTML report (``<out_dir>/<out_name>.json`` and ``.html``, plus
``<out_name>.tss_enrichment.csv`` when a TSS BED is given). The two plotted editing
distributions are also written as ``<out_name>.edits_per_fragment.csv`` and
``<out_name>.edit_rate_per_fragment.csv``, so they can be re-drawn without a rerun.
The HTML embeds the
multi-panel summary figure and documents the meaning of every metric inline.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import logging
import os
from collections import Counter, defaultdict
from collections.abc import Iterator
from datetime import datetime
from functools import partial
from typing import NamedTuple

import numpy as np
import pysam

from deamtools.utils import (
    Fragment,
    get_chrom_sizes_from_bam,
    get_version,
    iter_fragments,
    merge_fragment_bases,
    run_jobs,
)

logger = logging.getLogger(__name__)

# Histograms are stored as fixed-length arrays with a final overflow bin.
_MAX_EDITS = 100  # edits-per-fragment histogram: bins 0.._MAX_EDITS (last = overflow)
_MAX_FRAGLEN = 1000  # fragment-length histogram: bins 0.._MAX_FRAGLEN (overflow)
# Per-fragment edit-rate histogram: _RATE_BINS bins over [0, 1]. Real data piles
# up below 0.1, and the report plots this on a log-ish x axis, so 0.005-wide bins
# are needed to resolve the peak -- 0.02 bins leave it about five bars wide.
_RATE_BINS = 200
_MOTIF_WINDOW = 11  # bp window for the deaminase motif logo (odd; centre +/- 5)

# TSS-enrichment geometry, from the ENCODE ATAC-seq pipeline
# (encode_task_tss_enrich.py): a +/-2 kb window split into 400 bins, i.e. 10 bp
# each, with the outermost 100 bp on either side used as the background.
_TSS_BIN = 10
_TSS_EDGE = 100

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


# Fixed so that subsampled QC runs are reproducible without another CLI flag.
_SAMPLE_KEY = b"deamtools-qc-20260910"

# Library layout, decided from the head of the BAM before the main pass.
LAYOUT_PAIRED = "paired-end"
LAYOUT_SINGLE = "single-end"
_LAYOUT_PROBE = 10_000  # records read from the start of the file to decide it

# Companion CSVs: the numbers behind each plotted distribution, one file each,
# named <out_name>.<kind>.csv. The JSON names them rather than repeating them.
_CSV_EDITS = "edits_per_fragment"
_CSV_RATE = "edit_rate_per_fragment"
_CSV_TSS = "tss_enrichment"
_CSV_MOTIF = "motif_pfm"

# Y axis of the deaminase motif logo.
LOGO_BITS = "bits"
LOGO_FREQUENCY = "frequency"
LOGO_SCALES = (LOGO_BITS, LOGO_FREQUENCY)


def _csv_name(out_name: str, kind: str) -> str:
    return f"{out_name}.{kind}.csv"


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


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


def _zero_pair() -> list[int]:
    return [0, 0]


class _Stats:
    """Accumulator for QC counts; instances merge with :meth:`update`."""

    def __init__(self) -> None:
        self.qcfail = 0
        self.low_mapq = 0
        self.edit_counts: Counter[int] = Counter()
        self.length_counts: Counter[int] = Counter()
        self.direction_counts: Counter[tuple[int, int]] = Counter()
        self.zero_strata: Counter[tuple[str, bool]] = Counter()
        self.diagnostics: Counter[str] = Counter()
        self.background_pwm = np.zeros((_MOTIF_WINDOW, 4), dtype=np.int64)
        self.background_events = 0
        self.total = 0
        self.unmapped = 0
        self.duplicate = 0
        self.secondary = 0
        self.supplementary = 0
        self.proper_pair = 0
        self.passing = 0
        # Fragments, not records: a mate pair is merged into one (see
        # _accumulate_fragment), so these are the denominators for every
        # editing metric below.
        self.fragments = 0
        self.fragments_from_pairs = 0
        self.total_opportunities = 0
        self.total_edits = 0
        self.edits_per_fragment = np.zeros(_MAX_EDITS + 1, dtype=np.int64)
        self.fraglen = np.zeros(_MAX_FRAGLEN + 1, dtype=np.int64)
        # Per-fragment edit rate (edited C/G over editable C/G): distribution.
        self.edit_rate_hist = np.zeros(_RATE_BINS, dtype=np.int64)
        self.edit_rate_sum = 0.0
        self.edit_rate_n = 0
        # Deaminase motif: per-position A/C/G/T counts over the window around
        # each editing event (centre excluded, G-edits reverse-complemented).
        self.motif_pwm = np.zeros((_MOTIF_WINDOW, 4), dtype=np.int64)
        self.motif_events = 0
        # context -> [edits, opportunities]
        # A named factory, not a lambda: _Stats is pickled back from worker
        # processes, and a lambda default_factory cannot be pickled.
        self.context: dict[str, list[int]] = defaultdict(_zero_pair)

    def update(self, other: _Stats) -> None:
        peak = max(
            self.diagnostics["pending_mates_peak"],
            other.diagnostics["pending_mates_peak"],
        )
        self.qcfail += other.qcfail
        self.low_mapq += other.low_mapq
        for name in (
            "edit_counts",
            "length_counts",
            "direction_counts",
            "zero_strata",
            "diagnostics",
        ):
            getattr(self, name).update(getattr(other, name))
        self.diagnostics["pending_mates_peak"] = peak
        self.background_pwm += other.background_pwm
        self.background_events += other.background_events
        self.total += other.total
        self.unmapped += other.unmapped
        self.duplicate += other.duplicate
        self.secondary += other.secondary
        self.supplementary += other.supplementary
        self.proper_pair += other.proper_pair
        self.passing += other.passing
        self.fragments += other.fragments
        self.fragments_from_pairs += other.fragments_from_pairs
        self.total_opportunities += other.total_opportunities
        self.total_edits += other.total_edits
        self.edits_per_fragment += other.edits_per_fragment
        self.fraglen += other.fraglen
        self.edit_rate_hist += other.edit_rate_hist
        self.edit_rate_sum += other.edit_rate_sum
        self.edit_rate_n += other.edit_rate_n
        self.motif_pwm += other.motif_pwm
        self.motif_events += other.motif_events
        for ctx, (e, o) in other.context.items():
            slot = self.context[ctx]
            slot[0] += e
            slot[1] += o


def _keep_fragment(qname: str, fraction: float) -> bool:
    """Keep this fragment when subsampling? Decided from the read name.

    The decision has to be per *fragment*, not per record: drawing the two mates
    independently would leave most surviving fragments with only one mate, which
    would roughly halve every per-fragment edit count. Hashing the name gives
    both mates the same answer without having to remember any decisions.

    ``blake2b`` rather than the builtin ``hash()``, which is salted per process
    and would make reruns disagree; the seed goes in as the key.
    """
    digest = hashlib.blake2b(qname.encode(), digest_size=8, key=_SAMPLE_KEY).digest()
    return int.from_bytes(digest, "big") < fraction * (1 << 64)


def _accumulate_fragment(
    stats: _Stats,
    reads: Fragment,
    ref_seq: str,
    ref_len: int,
    min_baseq: int,
) -> None:
    """Fold one fragment (one record, or a merged mate pair) into ``stats``."""
    stats.fragments += 1
    if len(reads) > 1:
        stats.fragments_from_pairs += 1

    bases = merge_fragment_bases(reads, min_baseq, diagnostics=stats.diagnostics)
    ct = ga = 0

    frag_editable = 0  # distinct reference C/G covered (strand-agnostic)
    frag_edited = 0  # of those, how many show C->T or G->A
    for rpos, read_base in bases.items():
        ref_base = ref_seq[rpos]
        if read_base not in "ACGT" or ref_base not in "ACGT":
            continue
        stats.diagnostics["aligned_bases"] += 1
        if read_base != ref_base and (ref_base, read_base) not in (
            ("C", "T"),
            ("G", "A"),
        ):
            stats.diagnostics["non_edit_mismatches"] += 1
        ct += int(ref_base == "C" and read_base == "T")
        ga += int(ref_base == "G" and read_base == "A")

        # Per-fragment edit rate: every reference C or G is an editable base;
        # a C->T or G->A mismatch is an edit. Counted strand-agnostically and
        # without the flank requirement used for context below.
        if ref_base == "C":
            frag_editable += 1
            if read_base == "T":
                frag_edited += 1
        elif ref_base == "G":
            frag_editable += 1
            if read_base == "A":
                frag_edited += 1

        if rpos == 0 or rpos >= ref_len - 1:
            continue  # need flanking bases for the trinucleotide context

        # Strand-agnostic edit calling (matching bam2bw): a reference C may be
        # edited C->T and a reference G may be edited G->A, regardless of read
        # orientation. The G-centred context is reverse-complemented so both
        # are reported as C->T.
        if ref_base == "C":
            ctx = ref_seq[rpos - 1 : rpos + 2]
            is_edit = read_base == "T"
        elif ref_base == "G":
            ctx = _revcomp(ref_seq[rpos - 1 : rpos + 2])
            is_edit = read_base == "A"
        else:
            continue

        if any(b not in "ACGT" for b in ctx):
            continue

        # Context statistics have their own restricted denominator.
        slot = stats.context[ctx]
        slot[1] += 1
        if is_edit:
            slot[0] += 1

        # Opportunity-matched background and edited windows use identical bounds.
        half = _MOTIF_WINDOW // 2
        lo, hi = rpos - half, rpos + half + 1
        if lo >= 0 and hi <= ref_len:
            window = ref_seq[lo:hi]
            if ref_base == "G":
                window = _revcomp(window)
            if all(b in "ACGT" for b in window):
                stats.background_events += 1
                stats.motif_events += int(is_edit)
                for j, b in enumerate(window):
                    stats.background_pwm[j, _BASE_IDX[b]] += 1
                    if is_edit and j != half:
                        stats.motif_pwm[j, _BASE_IDX[b]] += 1

    frag_edits = frag_edited
    stats.total_opportunities += frag_editable
    stats.total_edits += frag_edited
    stats.edit_counts[frag_edits] += 1
    stats.direction_counts[(ct, ga)] += 1
    bounds = (0, 10, 25, 50, 100)
    labels = ("0", "1-10", "11-25", "26-50", "51-100", "101+")
    group = next(
        (labels[i] for i, n in enumerate(bounds) if frag_editable <= n), labels[-1]
    )
    stats.zero_strata[(group, frag_edits == 0)] += 1

    stats.edits_per_fragment[min(frag_edits, _MAX_EDITS)] += 1

    if frag_editable > 0:
        rate = frag_edited / frag_editable
        bin_idx = min(int(rate * _RATE_BINS), _RATE_BINS - 1)
        stats.edit_rate_hist[bin_idx] += 1
        stats.edit_rate_sum += rate
        stats.edit_rate_n += 1


def _detect_layout(bam_path: str, n: int = _LAYOUT_PROBE) -> str:
    """Is this a paired-end or a single-end library? Decided from ``n`` records.

    Paired-end if any of the first ``n`` records carries the paired flag
    (0x1). One is enough in both directions: a single-end library never sets
    it, and a paired-end library sets it on essentially every record --
    unmapped reads and orphans included, so it does not matter what sorts to
    the top of the file. The records are read from the start of the file
    rather than through the index, so unplaced reads are seen too.

    An empty BAM is reported as single-end, which only means that no
    pair-dependent metric is produced for it.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.head(n):
            if read.is_paired:
                return LAYOUT_PAIRED
    return LAYOUT_SINGLE


def _process_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    min_mapq: int,
    min_baseq: int,
    sample_fraction: float = 1.0,
) -> _Stats:
    """Accumulate read, editing, context and fragment-length stats for one chrom.

    Editing is accumulated per **fragment**: mates are held until their partner
    arrives and then merged, so a position both mates cover counts once. Records
    whose mate never arrives -- unpaired, mate filtered out, mate on another
    contig -- are folded in on their own at the end of the pass.

    When ``sample_fraction`` is below 1, a fragment is kept with that
    probability and skipped before any analysis, which is where the cost is.
    The draw is seeded and keyed on the read name, so a rerun samples the same
    fragments and both mates always share their fragment's fate.
    """
    stats = _Stats()
    subsampling = sample_fraction < 1.0

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        ref_seq = fasta.fetch(chrom).upper()
        ref_len = len(ref_seq)

        def passing_reads() -> Iterator[pysam.AlignedSegment]:
            """Record-level bookkeeping, yielding only what reaches a fragment."""
            for read in bam.fetch(chrom):
                if subsampling and not _keep_fragment(
                    read.query_name or "", sample_fraction
                ):
                    continue
                stats.total += 1
                if read.is_unmapped:
                    stats.unmapped += 1
                if read.is_duplicate:
                    stats.duplicate += 1
                if read.is_secondary:
                    stats.secondary += 1
                if read.is_supplementary:
                    stats.supplementary += 1
                stats.qcfail += int(read.is_qcfail)
                stats.low_mapq += int(
                    not read.is_unmapped and read.mapping_quality < min_mapq
                )

                if not _passes_filters(read, min_mapq):
                    continue
                stats.passing += 1
                for op, length in read.cigartuples or []:
                    if op in (0, 1, 4, 7, 8):
                        stats.diagnostics["query_bases"] += length
                    if op in (0, 2, 7, 8):
                        stats.diagnostics["reference_span_bases"] += length
                    if op in (1, 2, 4):
                        stats.diagnostics[
                            {
                                1: "inserted_bases",
                                2: "deleted_bases",
                                4: "soft_clipped_bases",
                            }[op]
                        ] += length

                # Fragment length from properly-paired read1 only, so a pair is
                # counted once.
                if read.is_proper_pair:
                    stats.proper_pair += 1
                    if read.is_read1 and read.template_length:
                        stats.length_counts[abs(read.template_length)] += 1
                        flen = min(abs(read.template_length), _MAX_FRAGLEN)
                        stats.fraglen[flen] += 1

                yield read

        for fragment in iter_fragments(
            passing_reads(), coordinate_sorted=True, diagnostics=stats.diagnostics
        ):
            _accumulate_fragment(stats, fragment, ref_seq, ref_len, min_baseq)

    return stats


class _TssResult(NamedTuple):
    """Aggregate TSS profile and its enrichment score.

    ``positions`` are bin centres in bp relative to the TSS, ``counts`` the raw
    insertion counts summed over all TSS, and ``profile`` the ENCODE-normalised
    signal (mean insertions per TSS divided by ``background``).
    """

    score: float
    positions: np.ndarray
    counts: np.ndarray
    profile: np.ndarray
    background: float
    n_tss: int
    bin_size: int


def _tss_sites(
    tss_path: str,
    chrom_sizes: dict[str, int],
    span: int,
    diagnostics: Counter[str] | None = None,
) -> Iterator[tuple[str, int, bool]]:
    """Yield ``(chrom, centre, is_minus)`` for TSS whose window fits the contig.

    The TSS is the midpoint of the BED record. Point annotations with
    start == end, 1-bp sites, and wider intervals are accepted. Strand is column 6 (BED6); records without it are treated as
    ``+``, except for the common four-column ``chrom start end strand`` TSS
    files, where column 4 is read as the strand instead.
    """
    with open(tss_path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            if diagnostics is not None:
                diagnostics["input_records"] += 1
            parts = line.rstrip("\n").split("\t")
            try:
                start, end = int(parts[1]), int(parts[2])
                if start < 0 or end < start:
                    raise ValueError("invalid interval")
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid TSS BED record: {line.strip()}") from exc
            chrom = parts[0]
            if chrom not in chrom_sizes:
                if diagnostics is not None:
                    diagnostics["unknown_contig"] += 1
                continue
            center = (start + end) // 2
            if center - span < 0 or center + span > chrom_sizes[chrom]:
                if diagnostics is not None:
                    diagnostics["outside_contig"] += 1
                continue
            if diagnostics is not None:
                diagnostics["used_records"] += 1
            strand_col = 3 if len(parts) == 4 else 5
            is_minus = len(parts) > strand_col and parts[strand_col].strip() == "-"
            yield chrom, center, is_minus


def _tss_enrichment(
    bam_path: str,
    tss_path: str,
    chrom_sizes: dict[str, int],
    min_mapq: int,
    flank: int,
    diagnostics: Counter[str] | None = None,
    sites: list[tuple[str, int, bool]] | None = None,
) -> _TssResult | None:
    """TSS enrichment, following the ENCODE ATAC-seq pipeline definition.

    Mirrors ``encode_task_tss_enrich.py`` of the ENCODE ATAC-seq pipeline:

    1. Take a ``+/-flank`` window around each TSS (ENCODE: 2 kb) and bin it at
       ``_TSS_BIN`` bp (ENCODE: 400 bins over 4 kb, i.e. 10 bp).
    2. Count insertion sites -- a read's 5' end, ``reference_start`` forward and
       ``reference_end - 1`` reverse -- into those bins, flipping the window for
       minus-strand TSS so upstream is always on the left.
    3. Average over TSS, then apply the Greenleaf normalisation: divide by the
       mean of the two edge means, each taken over the outermost ``_TSS_EDGE``
       bp (ENCODE: 100 bp), so the flanks sit at 1.
    4. The score is the **peak** of that normalised profile.

    The one deviation from ENCODE is deliberate. ENCODE reaches the insertion
    site indirectly, asking metaseq for read *coverage* shifted by
    ``-read_len/2`` so that each read's interval is centred on its cut site;
    that spreads every insertion over a read-length-wide box, which smooths the
    profile and depresses the peak. Counting the cut site itself at 1 bp is what
    that shift is trying to approximate, so scores here run slightly above
    ENCODE's for the same library, by more the longer the reads.

    Returns ``None`` when no TSS window fits the reference.
    """
    # A window narrower than one bin would leave nothing to bin; fall back to
    # per-bp resolution so small references (and tests) still produce a profile.
    bin_size = _TSS_BIN if flank >= _TSS_BIN else 1
    half_bins = flank // bin_size
    if half_bins == 0:
        logger.warning(f"  tss_flank={flank} is too small for a profile; skipping")
        return None
    n_bins = 2 * half_bins
    span = half_bins * bin_size  # effective flank after rounding to whole bins

    counts = np.zeros(n_bins, dtype=np.float64)
    n_tss = 0

    sites = sorted(
        _tss_sites(tss_path, chrom_sizes, span, diagnostics) if sites is None else sites
    )
    # Merge nearby windows for I/O only. Each original TSS (including duplicates)
    # still contributes its own profile. Cap a batch to avoid dense-contig growth.
    batches: list[list[tuple[str, int, bool]]] = []
    for site in sites:
        if (
            not batches
            or site[0] != batches[-1][-1][0]
            or site[1] - batches[-1][-1][1] > 2 * span
            or site[1] - batches[-1][0][1] > 1_000_000
        ):
            batches.append([])
        batches[-1].append(site)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for batch in batches:
            chrom = batch[0][0]
            start, end = batch[0][1] - span, batch[-1][1] + span
            cuts = []
            for read in bam.fetch(chrom, start, end):
                if not _passes_filters(read, min_mapq) or read.reference_end is None:
                    continue
                cut = (
                    read.reference_end - 1 if read.is_reverse else read.reference_start
                )
                if start <= cut < end:
                    cuts.append(cut)
            positions_in_batch = np.sort(np.asarray(cuts, dtype=np.int64))
            for _, center, is_minus in batch:
                lo, hi = np.searchsorted(
                    positions_in_batch, [center - span, center + span]
                )
                bins = (positions_in_batch[lo:hi] - center + span) // bin_size
                per_tss = np.bincount(bins, minlength=n_bins)
                counts += per_tss[::-1] if is_minus else per_tss
                n_tss += 1

    positions = (np.arange(n_bins) - half_bins) * bin_size + bin_size / 2.0

    if n_tss == 0:
        logger.warning("  no usable TSS found; skipping TSS enrichment")
        return None

    mean_per_tss = counts / n_tss
    edge_bins = max(1, min(_TSS_EDGE // bin_size, half_bins))
    # ENCODE averages the two edge means rather than pooling their bins. With
    # equal-width edges the two agree; keep the reference form regardless.
    background = float(
        (mean_per_tss[:edge_bins].mean() + mean_per_tss[-edge_bins:].mean()) / 2.0
    )
    if background <= 0:
        logger.warning(
            "  no insertions in the TSS flank background; enrichment undefined"
        )
        return _TssResult(
            float("nan"),
            positions,
            counts,
            np.full(n_bins, np.nan),
            0.0,
            n_tss,
            bin_size,
        )

    profile = mean_per_tss / background
    score = float(profile.max())
    logger.info(
        f"  TSS enrichment = {score:.2f} (peak of the normalised profile over "
        f"{n_tss} TSS, {bin_size}-bp bins, +/-{span} bp)"
    )
    return _TssResult(score, positions, counts, profile, background, n_tss, bin_size)


def _write_edits_per_fragment_csv(path: str, hist: np.ndarray) -> None:
    """Write the edits-per-fragment histogram, the numbers behind panel 2.

    One row per edit count ``0 .. _MAX_EDITS``. The last row is an overflow bin
    holding every fragment with *at least* that many edits, and is flagged in
    ``is_overflow`` so a re-plot can label or drop it.
    """
    total = int(hist.sum())
    last = len(hist) - 1
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["edits", "fragments", "fraction", "is_overflow"])
        for edits, count in enumerate(hist):
            fraction = count / total if total else 0.0
            writer.writerow([edits, int(count), f"{fraction:.6g}", edits == last])


def _write_edit_rate_csv(path: str, hist: np.ndarray) -> None:
    """Write the per-fragment edit-rate histogram, the numbers behind panel 3.

    One row per bin of width ``1 / _RATE_BINS``. Bins are half-open
    ``[bin_start, bin_end)`` except the last, which also holds a rate of
    exactly 1.
    """
    n_bins = len(hist)
    total = int(hist.sum())
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bin_start", "bin_end", "fragments", "fraction"])
        for i, count in enumerate(hist):
            fraction = count / total if total else 0.0
            writer.writerow(
                [
                    f"{i / n_bins:g}",
                    f"{(i + 1) / n_bins:g}",
                    int(count),
                    f"{fraction:.6g}",
                ]
            )


def _write_tss_csv(path: str, tss: _TssResult) -> None:
    """Write the aggregate TSS profile -- the numbers behind the plot -- as CSV."""
    mean_per_tss = tss.counts / tss.n_tss
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["position", "insertions", "mean_insertions_per_tss", "normalized"]
        )
        for pos, count, mean, norm in zip(
            tss.positions, tss.counts, mean_per_tss, tss.profile, strict=True
        ):
            writer.writerow(
                [
                    f"{pos:g}",
                    int(count),
                    f"{mean:.6g}",
                    f"{float(norm):.6g}" if np.isfinite(norm) else "",
                ]
            )


def _count_summary(counts: dict[int, int]) -> dict:
    """Exact mean and conventional median from an untruncated frequency table."""
    n = sum(counts.values())
    if not n:
        return {"n": 0, "mean": None, "median": None}
    ranks = ((n - 1) // 2, n // 2)
    cumulative = 0
    middle: list[int] = []
    for value, count in sorted(counts.items()):
        middle.extend(
            value for rank in ranks if cumulative <= rank < cumulative + count
        )
        cumulative += count
    return {
        "n": n,
        "mean": sum(k * v for k, v in counts.items()) / n,
        "median": sum(middle) / 2,
    }


def _histogram_summary(hist: np.ndarray) -> dict:
    return _count_summary({i: int(n) for i, n in enumerate(hist) if n})


def _ratio(a: int | float, b: int | float) -> float | None:
    return a / b if b else None


def _build_metrics(
    stats: _Stats,
    tss: _TssResult | None,
    out_name: str,
    layout: str = LAYOUT_PAIRED,
) -> dict:
    edit_rate = (
        stats.total_edits / stats.total_opportunities
        if stats.total_opportunities
        else None
    )
    per_fragment = _count_summary(stats.edit_counts)
    fraglen = _count_summary(stats.length_counts)

    # Per-fragment edit rate (edited C/G over editable C/G). Mean is exact;
    # median is taken from the histogram bin centres.
    rate_mean = stats.edit_rate_sum / stats.edit_rate_n if stats.edit_rate_n else None
    if stats.edit_rate_n:
        bin_median = _histogram_summary(stats.edit_rate_hist)["median"]
        rate_median = (bin_median + 0.5) / _RATE_BINS
    else:
        rate_median = None

    context = {
        ctx: {
            "edits": e,
            "opportunities": o,
            "edit_fraction": (e / o if o else 0.0),
        }
        for ctx, (e, o) in sorted(stats.context.items())
    }

    metrics: dict = {
        "schema_version": "2.0",
        "library_layout": layout,
        "context_summary": {
            "total_edits": sum(e for e, _ in stats.context.values()),
            "total_opportunities": sum(o for _, o in stats.context.values()),
            "definition": "ACGT trinucleotide windows; separate from unrestricted global counts",
        },
        "diagnostics": dict(stats.diagnostics),
        "reads": {
            "scope": "sampled records, including unplaced unmapped reads",
            "filter_counts_overlap": True,
            "qcfail": stats.qcfail,
            "low_mapq": stats.low_mapq,
            "total": stats.total,
            "passing": stats.passing,
            "unmapped": stats.unmapped,
            "duplicate": stats.duplicate,
            "duplicate_rate": _ratio(stats.duplicate, stats.total),
            "secondary": stats.secondary,
            "supplementary": stats.supplementary,
        },
        "fragments": {
            "total": stats.fragments,
            "from_mate_pairs": stats.fragments_from_pairs,
            "from_single_records": stats.fragments - stats.fragments_from_pairs,
        },
        "editing": {
            "total_opportunities": stats.total_opportunities,
            "total_edits": stats.total_edits,
            "global_edit_rate": edit_rate,
            "status": "ok" if stats.total_opportunities else "no_editable_bases",
            "zero_edit_fragments": stats.edit_counts[0],
            "zero_edit_fraction": _ratio(stats.edit_counts[0], stats.fragments),
            "zero_edit_by_opportunities_csv": _csv_name(
                out_name, "zero_edit_by_opportunities"
            ),
            "direction_joint_csv": _csv_name(out_name, "edit_directions"),
            "ct_edits": sum(ct * n for (ct, ga), n in stats.direction_counts.items()),
            "ga_edits": sum(ga * n for (ct, ga), n in stats.direction_counts.items()),
            "both_direction_fragments": sum(
                n for (ct, ga), n in stats.direction_counts.items() if ct and ga
            ),
            "mean_edits_per_fragment": per_fragment["mean"],
            "median_edits_per_fragment": per_fragment["median"],
            "edits_per_fragment_csv": _csv_name(out_name, _CSV_EDITS),
        },
        "edit_rate_per_fragment": {
            # Not called `n_reads`: that is the --n_reads subsample size, which
            # is a plain uniform draw with no condition on the read at all.
            "n_fragments_with_editable_bases": stats.edit_rate_n,
            "mean": rate_mean,
            "median": rate_median,
            "median_method": "approximate, 0.005-wide bin centres",
            # Like the TSS profile, the histogram's one home is its CSV.
            "histogram_csv": _csv_name(out_name, _CSV_RATE),
        },
        "context": context,
        "motif": {
            "window": _MOTIF_WINDOW,
            "n_events": stats.motif_events,
            "pfm_csv": _csv_name(out_name, _CSV_MOTIF),
            "background_events": stats.background_events,
            "enrichment_csv": _csv_name(out_name, "motif_enrichment"),
        },
    }
    # Pair-dependent metrics exist only for a paired-end library. For a
    # single-end one they are left out rather than reported as 0: a
    # proper-pair rate of 0 reads as a mapping failure, and a mean fragment
    # length of 0 bp is not a length at all -- there is no insert to measure.
    if layout == LAYOUT_PAIRED:
        metrics["reads"]["proper_pair"] = stats.proper_pair
        metrics["reads"]["proper_pair_rate"] = _ratio(stats.proper_pair, stats.passing)
        metrics["fragment_length"] = {
            "mean": fraglen["mean"],
            "median": fraglen["median"],
            "n_pairs": fraglen["n"],
            "histogram_csv": _csv_name(out_name, "fragment_length"),
        }
    if tss is not None:
        # The profile itself is not duplicated here: it is a few hundred numbers
        # and its one home is the CSV, named below so the JSON still points at it.
        metrics["tss_enrichment"] = {
            "score": tss.score if np.isfinite(tss.score) else None,
            "status": "ok" if np.isfinite(tss.score) else "zero_background",
            "n_tss": tss.n_tss,
            "flank": int(len(tss.counts) // 2 * tss.bin_size),
            "bin_size": tss.bin_size,
            "background": tss.background,
            "total_insertions": int(tss.counts.sum()),
            "profile_csv": _csv_name(out_name, _CSV_TSS),
        }
    for numerator, denominator, label in (
        ("non_edit_mismatches", "aligned_bases", "non_edit_mismatch_rate"),
        ("overlap_conflicts", "overlap_positions", "overlap_conflict_rate"),
        ("equal_quality_conflicts", "overlap_positions", "equal_quality_conflict_rate"),
        ("inserted_bases", "query_bases", "insertion_rate"),
        ("soft_clipped_bases", "query_bases", "soft_clip_rate"),
        ("deleted_bases", "reference_span_bases", "deletion_rate"),
    ):
        metrics["diagnostics"][numerator] = stats.diagnostics[numerator]
        metrics["diagnostics"][denominator] = stats.diagnostics[denominator]
        metrics["diagnostics"][label] = _ratio(
            stats.diagnostics[numerator], stats.diagnostics[denominator]
        )
    return metrics


def _motif_counts_df(pwm: np.ndarray, n_events: int):
    """The motif as a position x base count matrix, the target cytosine included.

    ``pwm`` is accumulated with the edited base itself left out, so its centre
    row is all zeros. That position is known, though: every event is centred on
    a reference C once G->A events are reverse-complemented into the C->T
    orientation, so the centre row is filled in as ``C = n_events``. Flank rows
    sum to at most ``n_events`` -- a flank base that is N is not counted.

    Returns a DataFrame indexed by offset from the edited base (``position``,
    centre 0) with integer columns ``A C G T``: what ``<out_name>.motif_pfm.csv``
    holds, and what the frequency logo normalises.
    """
    import pandas as pd

    window = pwm.shape[0]
    counts = pd.DataFrame(pwm.astype(np.int64), columns=["A", "C", "G", "T"])
    counts.index = pd.Index(np.arange(window) - window // 2, name="position")
    counts.loc[0] = [0, n_events, 0, 0]
    return counts


def _write_motif_csv(path: str, pwm: np.ndarray, n_events: int) -> None:
    """Write the motif count matrix, the numbers behind the logo."""
    _motif_counts_df(pwm, n_events).to_csv(path)


def _motif_bits_df(pwm: np.ndarray):
    """Convert a per-position A/C/G/T count PWM to information content (bits).

    Returns a DataFrame indexed by position offset (centre = 0) with one column
    per base, ready for :class:`logomaker.Logo`.
    """
    import pandas as pd

    window = pwm.shape[0]
    df = pd.DataFrame(pwm, columns=["A", "C", "G", "T"], dtype=float)
    totals = df.sum(axis=1).replace(0, np.nan)
    probs = df.div(totals, axis=0).fillna(0.0)
    p = probs.to_numpy()
    log_p = np.zeros_like(p)
    with np.errstate(divide="ignore", invalid="ignore"):
        np.log2(p, where=(p > 0), out=log_p)
    h = -np.sum(p * log_p, axis=1)
    info = np.log2(p.shape[1]) - h  # log2(4) - entropy
    bits = probs.multiply(info, axis=0)
    bits.index = np.arange(window) - window // 2
    bits.index.name = "position"
    return bits


def _motif_logo_base64(
    pwm: np.ndarray, n_events: int = 0, scale: str = LOGO_FREQUENCY
) -> str | None:
    """Render the deaminase motif logo; base64 PNG, or None if there are no events.

    ``scale="bits"`` plots information content, with the edited base left out:
    it is always C, so it would carry the full 2 bits and flatten the flanks,
    which rarely reach 0.15. ``scale="frequency"`` plots each base's frequency
    at each offset on a 0-1 axis, with the target C drawn at position 0.
    """
    if pwm.sum() == 0:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import logomaker  # noqa: E402
    import matplotlib.pyplot as plt  # noqa: E402

    if scale == LOGO_FREQUENCY:
        counts = _motif_counts_df(pwm, n_events)
        freq = counts.div(counts.sum(axis=1), axis=0).fillna(0.0)  # count -> freq
        fig, ax = plt.subplots(figsize=(4, 2.5))
        logo = logomaker.Logo(
            freq,
            ax=ax,
            color_scheme="classic",  # A green, C blue, G orange, T red
            show_spines=False,
        )
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Distance from target cytosine")
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        logo.style_spines(visible=False)
        logo.style_spines(spines=["left", "bottom"], visible=True, linewidth=0.8)
        ax.tick_params(labelsize=7)
    else:
        bits = _motif_bits_df(pwm)
        fig, ax = plt.subplots(figsize=(7, 2.4))
        logo = logomaker.Logo(bits, ax=ax, baseline_width=0)
        logo.style_spines(visible=False)
        logo.style_spines(spines=["left", "bottom"], visible=True)
        ax.set_xlabel("distance from edited base")
        ax.set_ylabel("bits")
        ax.set_title("Deaminase sequence motif")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


# Per-metric documentation shown in the HTML report. Each entry maps a JSON
# field to a plain-language explanation of what it means and how to read it.
_METRIC_DOCS: dict[str, dict[str, str]] = {
    "reads": {
        "scope": "Counters describe sampled records; file_reads is unsampled.",
        "filter_counts_overlap": "Flag and MAPQ reasons overlap; do not sum them as disjoint failures.",
        "qcfail": "Records carrying the QC-fail flag (0x200).",
        "low_mapq": "Mapped records below min_mapq; may also fail other filters.",
        "total": "Sampled read records before filtering, including unplaced unmapped records. See file_reads for full-file counts.",
        "passing": (
            "Reads kept after filtering (primary, non-duplicate, not QC-fail, "
            "MAPQ &ge; min_mapq). All editing, context and fragment metrics are "
            "computed from these reads only."
        ),
        "unmapped": "Reads flagged unmapped (SAM flag 0x4).",
        "duplicate": "Reads flagged as PCR/optical duplicates (flag 0x400).",
        "duplicate_rate": (
            "Duplicate-flagged records as a fraction of analyzed records. "
            "This may reflect upstream deduplication; it does not establish original library complexity."
        ),
        "secondary": "Secondary alignments (flag 0x100).",
        "supplementary": "Supplementary / chimeric alignments (flag 0x800).",
        "proper_pair": "Passing reads flagged as a proper pair (flag 0x2).",
        "proper_pair_rate": (
            "Proper pairs as a fraction of passing reads. Low values can signal "
            "insert-size or mapping problems."
        ),
    },
    "fragments": {
        "total": (
            "Fragments the editing metrics were computed over. Every metric "
            "below is per fragment, not per read: a mate pair is merged first, "
            "so a reference position both mates cover counts once."
        ),
        "from_mate_pairs": "Fragments built by merging two mates.",
        "from_single_records": (
            "Fragments that were a single record &mdash; unpaired reads, and "
            "reads whose mate was unmapped, filtered out, or on another contig."
        ),
    },
    "editing": {
        "status": "ok or no_editable_bases; undefined rates are n/a.",
        "zero_edit_fragments": "Fragments with zero edits, including fragments with no editable bases.",
        "zero_edit_fraction": "Zero-edit fragments / all analyzed fragments; not a closed-chromatin call.",
        "ct_edits": "C-to-T events independent of read orientation.",
        "ga_edits": "G-to-A events independent of read orientation.",
        "both_direction_fragments": "Fragments carrying at least one event of each direction.",
        "direction_joint_csv": "Joint C-to-T/G-to-A edit counts per fragment, as a sparse frequency table.",
        "zero_edit_by_opportunities_csv": "Zero-edit fraction stratified by covered editable bases.",
        "total_opportunities": (
            "Editable reference positions covered by passing fragments: any "
            "reference C or G (counted regardless of read orientation), with "
            "an A/C/G/T read base and base quality &ge; min_baseq. A "
            "position covered by both mates counts once."
        ),
        "total_edits": (
            "Opportunities showing a deamination event &mdash; a C&rarr;T "
            "mismatch at a reference C or a G&rarr;A mismatch at a reference G, "
            "regardless of read orientation (matching bam2bw)."
        ),
        "global_edit_rate": (
            "Overall fraction of covered editable C/G positions showing C→T or G→A editing. "
            "Interpret relative to matched controls and library conditions."
        ),
        "mean_edits_per_fragment": (
            "Average number of edits per fragment. Deaminase fragments carry "
            "many edits, unlike the two Tn5 cut sites of a standard ATAC read."
        ),
        "median_edits_per_fragment": "Median number of edits per fragment.",
        "edits_per_fragment_csv": (
            "File holding the edits-per-fragment histogram plotted above (one "
            "row per edit count; the last row is an overflow bin)."
        ),
    },
    "edit_rate_per_fragment": {
        "n_fragments_with_editable_bases": (
            "Fragments covering at least one reference C or G, which is what a "
            "per-fragment edit rate needs a denominator for. Note this counts "
            "editable bases, not editing events: a fragment with no edit at all "
            "still contributes, at rate 0. Unrelated to <code>--n_reads</code>."
        ),
        "mean": (
            "Mean of the per-fragment edit rate, where a fragment's rate = "
            "(edited C/G) / (editable C/G), counted strand-agnostically over "
            "every distinct reference C and G the fragment covers. Normalising "
            "by the number of editable bases makes fragments comparable "
            "regardless of length or composition."
        ),
        "median": "Median of the per-fragment edit-rate distribution.",
        "histogram_csv": (
            "File holding the edit-rate histogram plotted above (one row per "
            "bin, with its start, end, fragment count and fraction)."
        ),
    },
    "fragment_length": {
        "mean": (
            "Mean fragment length from properly-paired read 1 (absolute template "
            "length). An ATAC-style library shows nucleosome laddering."
        ),
        "median": "Median fragment length.",
        "n_pairs": "Number of read pairs contributing a fragment length.",
    },
    "tss_enrichment": {
        "score": (
            "Peak of the flank-normalised insertion profile around TSS, as "
            "computed using unshifted 5-prime sites and 10-bp bins. Compare matched "
            "protocols and parameters; no universal pass/fail threshold is applied."
        ),
        "n_tss": (
            "Number of TSS that contributed, i.e. those whose full window fits "
            "inside the contig and whose contig is present in the BAM header."
        ),
        "flank": "Half-width in bp of the window around each TSS (ENCODE: 2000).",
        "bin_size": "Width in bp of one profile bin (ENCODE: 10).",
        "background": (
            "Mean insertions per TSS per bin in the outermost 100 bp on each "
            "side; the profile is divided by this, so the flanks sit at 1 "
            "(the Greenleaf normalisation ENCODE applies)."
        ),
        "total_insertions": "Insertion sites counted across all TSS windows.",
        "profile_csv": "File holding the per-bin numbers plotted above.",
    },
}


def _write_diagnostics(
    out_dir: str, out_name: str, stats: _Stats, metrics: dict
) -> None:
    """Write reusable diagnostics for reports and downstream analyses."""

    def write(kind, header, rows):
        with open(os.path.join(out_dir, _csv_name(out_name, kind)), "w") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)

    write(
        "edit_directions",
        ["ct_edits", "ga_edits", "fragments"],
        [(ct, ga, n) for (ct, ga), n in sorted(stats.direction_counts.items())],
    )
    strata = []
    for group in ("0", "1-10", "11-25", "26-50", "51-100", "101+"):
        zero = stats.zero_strata[(group, True)]
        total = zero + stats.zero_strata[(group, False)]
        strata.append((group, total, zero, _ratio(zero, total)))
    write(
        "zero_edit_by_opportunities",
        ["editable_bases", "fragments", "zero_edit_fragments", "zero_edit_fraction"],
        strata,
    )
    write("fragment_length", ["length", "pairs"], sorted(stats.length_counts.items()))
    write(
        "context",
        ["context", "edits", "opportunities", "edit_fraction"],
        [
            (ctx, d["edits"], d["opportunities"], d["edit_fraction"])
            for ctx, d in metrics["context"].items()
        ],
    )
    event_counts = stats.motif_pwm.copy()
    event_counts[_MOTIF_WINDOW // 2, _BASE_IDX["C"]] = stats.motif_events
    rows = []
    for j in range(_MOTIF_WINDOW):
        for base, bi in _BASE_IDX.items():
            event_freq = _ratio(int(event_counts[j, bi]), stats.motif_events)
            background_freq = _ratio(
                int(stats.background_pwm[j, bi]), stats.background_events
            )
            # No pseudocount: absent events/background are explicitly undefined.
            fold = (
                event_freq / background_freq
                if event_freq is not None and background_freq
                else None
            )
            rows.append(
                (
                    j - _MOTIF_WINDOW // 2,
                    base,
                    int(event_counts[j, bi]),
                    int(stats.background_pwm[j, bi]),
                    event_freq,
                    background_freq,
                    fold,
                )
            )
    write(
        "motif_enrichment",
        [
            "position",
            "base",
            "edited_count",
            "opportunity_count",
            "edited_frequency",
            "background_frequency",
            "fold_enrichment",
        ],
        rows,
    )


def run_qc(
    bam_path: str,
    fasta_path: str,
    out_dir: str,
    out_name: str,
    tss_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    threads: int = 1,
    tss_flank: int = 2000,
    plot: bool = True,
    n_reads: int | None = None,
    logo_scale: str = LOGO_FREQUENCY,
    report_dir: str | None = None,
    command: str | None = None,
    redact_paths: bool = False,
) -> dict:
    """Compute QC metrics for a deaminase chromatin-accessibility BAM.

    Writes a machine-readable ``<out_dir>/<out_name>.json`` and a
    self-contained, MultiQC-style ``<out_dir>/<out_name>.html`` report that
    contains independent plots, concise metric tables and QC interpretation. With
    ``tss_path``, the aggregate TSS profile behind the report's plot is also
    written to ``<out_dir>/<out_name>.tss_enrichment.csv``. The edits-per-fragment
    and edit-rate histograms behind the summary figure always go to
    ``<out_name>.edits_per_fragment.csv`` and ``<out_name>.edit_rate_per_fragment.csv``,
    written even when ``plot`` is False; the JSON names each file rather than
    repeating the numbers.

    Parameters
    ----------
    bam_path : str
        Coordinate-sorted, indexed BAM file.
    fasta_path : str
        Reference FASTA indexed with ``samtools faidx``.
    out_dir : str
        Output directory. Created if it does not exist.
    out_name : str
        Base name (without extension) for the ``.json`` and ``.html`` outputs.
    tss_path : str, optional
        BED file of transcription start sites. When supplied, an ENCODE-style
        TSS enrichment score and aggregate profile are computed. Column 6 is
        read as the strand when present, so minus-strand TSS are flipped.
    min_mapq : int, default 20
        Minimum read mapping quality.
    min_baseq : int, default 20
        Minimum base quality for a position to count as an editing opportunity.
    threads : int, default 1
        Number of worker processes; chromosomes are processed in parallel.
        With 1 everything runs in this process.
    tss_flank : int, default 2000
        Half-width (bp) of the window around each TSS for enrichment. The
        ENCODE default is 2000; it is rounded down to a whole number of 10-bp
        bins.
    plot : bool, default True
        Whether to render and embed the summary figure in the HTML report.
    logo_scale : {"frequency", "bits"}, default "frequency"
        Y axis of the deaminase motif logo. ``"frequency"`` plots each base's
        frequency per offset on a 0-1 axis with the target C at position 0.
        ``"bits"`` plots information content with the edited base left out (it
        is always C, and would take the full 2 bits); it is the more standard
        logo, but a deaminase's preference is weak -- the flanks stay under
        ~0.15 bits -- so the letters come out barely legible, which is why
        frequency is the default. The counts behind either are written to
        ``<out_name>.motif_pfm.csv`` regardless.
    n_reads : int, optional
        Subsample to approximately this many reads instead of using all of
        them, for a faster pass over a large BAM. Fragments are drawn uniformly
        across the genome (each kept with probability ``n_reads / total``,
        which is read from the BAM index, so nothing is scanned twice), and
        the draw is seeded, so a rerun samples the same fragments. Ignored when
        the BAM already holds fewer reads than this.

        The unit of the draw is the **fragment**, keyed on the read name, so
        both mates always share their fragment's fate. Drawing the two mates
        independently would leave most surviving fragments with only one mate
        and roughly halve every per-fragment edit count; keying on the name
        costs one hash per record and avoids that entirely.

        The draw is blind to what a read contains: the coin is flipped before
        the read is looked at, so reads carrying no editing event are kept at
        exactly the same rate as reads full of them. (Do not confuse this with
        ``edit_rate_per_fragment.n_fragments_with_editable_bases``, which *is*
        conditional -- on covering a reference C or G, not on being edited.)

        Rates and distributions -- editing rate, duplicate rate, context bias,
        the motif PWM, fragment lengths -- are unbiased under this sampling.
        The absolute counts in the report are counts *of the sample*, not
        estimates of the whole file, and the ``sampling`` block of the metrics
        records the fraction so they can be scaled if needed. TSS enrichment
        always uses every read in its windows, since it is a ratio computed
        over a small part of the genome and subsampling it would only add
        noise.

    report_dir : str, optional
        Write report.html, assets and CSVs to a portable report directory.
        The default is a single standalone HTML in out_dir.
    command : str, optional
        Recorded invocation, only when supplied by the caller.

    redact_paths : bool, default False
        Remove directory components from HTML and portable report provenance.
        Local QC JSON retains the full recorded provenance.

    Returns
    -------
    dict
        The metrics dictionary (also written to ``<out_dir>/<out_name>.json``).
    """
    if min_mapq < 0 or min_baseq < 0 or threads < 1 or tss_flank < 1:
        raise ValueError(
            "Quality thresholds must be nonnegative; threads and tss_flank positive"
        )
    if logo_scale not in LOGO_SCALES:
        raise ValueError(
            f"logo_scale must be one of {', '.join(LOGO_SCALES)}, got {logo_scale!r}"
        )
    # Fail on an invalid output location before scanning a large BAM.
    os.makedirs(out_dir, exist_ok=True)
    if report_dir is not None:
        os.makedirs(os.path.join(report_dir, "assets"), exist_ok=True)
    logger.info("Running qc")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_sizes = get_chrom_sizes_from_bam(bam)
        # The index carries per-contig counts, so the total costs no scan.
        total_reads = (
            sum(st.mapped + st.unmapped for st in bam.get_index_statistics())
            + bam.nocoordinate
        )
        file_counts = {
            "total": total_reads,
            "mapped": bam.mapped,
            "unmapped": bam.unmapped,
            "unplaced": bam.nocoordinate,
        }
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom, length in chrom_sizes.items():
            if (
                chrom not in fasta.references
                or fasta.get_reference_length(chrom) != length
            ):
                raise ValueError(f"BAM/FASTA contig mismatch: {chrom}")
    # Parse once before the expensive BAM scan, then reuse the validated sites.
    tss_diagnostics: Counter[str] = Counter()
    tss_sites = None
    if tss_path is not None:
        bin_size = _TSS_BIN if tss_flank >= _TSS_BIN else 1
        span = (tss_flank // bin_size) * bin_size
        tss_sites = list(_tss_sites(tss_path, chrom_sizes, span, tss_diagnostics))
        logger.info(f"Validated {len(tss_sites)} usable TSS before BAM processing")
    chroms = list(chrom_sizes.keys())

    layout = _detect_layout(bam_path)
    logger.info(
        f"Library layout: {layout} (from the first {_LAYOUT_PROBE:,} record(s))"
    )

    sample_fraction = 1.0
    if n_reads is not None:
        if n_reads <= 0:
            raise ValueError(f"n_reads must be positive, got {n_reads}")
        if total_reads and n_reads < total_reads:
            sample_fraction = n_reads / total_reads
            logger.info(
                f"Subsampling to ~{n_reads} of {total_reads} read(s) "
                f"(fraction {sample_fraction:.4g})"
            )
        else:
            logger.info(
                f"n_reads={n_reads} >= {total_reads} read(s) in the BAM; "
                "using all reads"
            )

    logger.info(f"Processing {len(chroms)} chromosome(s) with {threads} worker(s)")

    jobs = [
        partial(
            _process_chrom,
            bam_path,
            fasta_path,
            c,
            min_mapq,
            min_baseq,
            sample_fraction,
        )
        for c in chroms
    ]
    merged = _Stats()
    for chrom_stats in run_jobs(jobs, threads):
        merged.update(chrom_stats)

    # Unplaced records are absent from every chromosome fetch.
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch("*"):
            if sample_fraction < 1 and not _keep_fragment(
                read.query_name or "", sample_fraction
            ):
                continue
            merged.total += 1
            merged.unmapped += int(read.is_unmapped)
            merged.duplicate += int(read.is_duplicate)
            merged.secondary += int(read.is_secondary)
            merged.supplementary += int(read.is_supplementary)
            merged.qcfail += int(read.is_qcfail)

    logger.info(
        f"  {merged.passing} passing read(s) in {merged.fragments} fragment(s) "
        f"({merged.fragments_from_pairs} from merged mate pairs); "
        f"{merged.total_edits}/{merged.total_opportunities} edits/opportunities"
    )

    tss: _TssResult | None = None
    if tss_path is not None:
        logger.info(f"TSS:   {tss_path}")
        tss = _tss_enrichment(
            bam_path, tss_path, chrom_sizes, min_mapq, tss_flank, sites=tss_sites
        )

    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f"{out_name}.json")
    metrics = _build_metrics(merged, tss, out_name, layout)
    metrics["provenance"] = {
        "sample": out_name,
        "deamtools_version": get_version(),
        "generated_at": datetime.now().astimezone().isoformat(),
        "bam": os.path.abspath(bam_path),
        "fasta": os.path.abspath(fasta_path),
        "tss": os.path.abspath(tss_path) if tss_path else None,
        "parameters": {
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "threads": threads,
            "tss_flank": tss_flank,
            "n_reads": n_reads,
            "plot": plot,
            "logo_scale": logo_scale,
            "report_dir": report_dir,
            "redact_paths": redact_paths,
        },
    }
    if command is not None:
        metrics["provenance"]["command"] = command
    metrics["file_reads"] = file_counts
    metrics["status"] = {
        "editing": metrics["editing"]["status"],
        "tss": (
            "not_requested"
            if tss_path is None
            else (
                "no_usable_tss"
                if tss is None
                else "ok" if np.isfinite(tss.score) else "zero_background"
            )
        ),
        "tss_annotation": dict(tss_diagnostics),
    }
    metrics["sampling"] = {
        "subsampled": sample_fraction < 1.0,
        "requested_reads": n_reads,
        "fraction": sample_fraction,
        "reads_in_bam": total_reads,
    }
    metrics["motif"]["logo_scale"] = logo_scale
    metrics["context_summary"]["csv"] = _csv_name(out_name, "context")

    # The numbers behind each plot, written whether or not the plots are.
    edits_csv = os.path.join(out_dir, _csv_name(out_name, _CSV_EDITS))
    logger.info(f"Writing {edits_csv}")
    _write_edits_per_fragment_csv(edits_csv, merged.edits_per_fragment)

    rate_csv = os.path.join(out_dir, _csv_name(out_name, _CSV_RATE))
    logger.info(f"Writing {rate_csv}")
    _write_edit_rate_csv(rate_csv, merged.edit_rate_hist)

    if tss is not None:
        tss_csv = os.path.join(out_dir, _csv_name(out_name, _CSV_TSS))
        logger.info(f"Writing {tss_csv}")
        _write_tss_csv(tss_csv, tss)

    motif_csv = os.path.join(out_dir, _csv_name(out_name, _CSV_MOTIF))
    logger.info(f"Writing {motif_csv}")
    _write_motif_csv(motif_csv, merged.motif_pwm, merged.motif_events)

    logger.info(f"Writing {json_path}")
    with open(json_path, "w") as f:
        json.dump(metrics, f, indent=2, allow_nan=False)

    _write_diagnostics(out_dir, out_name, merged, metrics)
    from deamtools.qc.report import write_report

    motif_b64 = (
        _motif_logo_base64(merged.motif_pwm, merged.motif_events, logo_scale)
        if plot
        else None
    )
    html_path = write_report(
        metrics,
        out_dir,
        out_name,
        plot=plot,
        motif_b64=motif_b64,
        report_dir=report_dir,
        descriptions=_METRIC_DOCS,
        redact_paths=redact_paths,
    )
    logger.info(f"Wrote {html_path}")

    logger.info("Done")
    return metrics
