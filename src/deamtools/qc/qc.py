"""Quality-control metrics for deaminase-based chromatin accessibility data.

Summarises a coordinate-sorted BAM together with its reference FASTA into the
metrics most useful for judging a deaminase footprinting experiment:

* **Read statistics** — totals plus the fraction of duplicate, properly-paired,
  secondary and supplementary reads.
* **Editing statistics** — the genome-wide deamination rate (edits divided by
  the number of editable C/G *opportunities* covered by passing reads) and the
  distribution of edits per read. A high, accessibility-driven edit rate is the
  primary signal that the deaminase treatment worked.
* **Trinucleotide context bias** — the edit fraction broken down by the
  trinucleotide centred on the edited cytosine. Edits are called
  strand-agnostically (matching :mod:`deamtools.preprocessing.bam2bw`): any
  reference ``C->T`` or ``G->A`` mismatch counts, and ``G``-centred contexts
  are reverse-complemented so both are reported in the unified ``C``-centred
  orientation. This is
  the enzyme's sequence-preference fingerprint (e.g. DddA's ``TC`` preference).
* **Fragment-length distribution** — from the template length of properly-paired
  read pairs.
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
``<out_name>.tss_enrichment.csv`` when a TSS BED is given). The HTML embeds the
multi-panel summary figure and documents the meaning of every metric inline.
"""

from __future__ import annotations

import base64
import csv
import io
import json
import logging
import os
import random
from collections import defaultdict
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import NamedTuple

import numpy as np
import pysam

from deamtools.utils import get_chrom_sizes_from_bam

logger = logging.getLogger(__name__)

# Histograms are stored as fixed-length arrays with a final overflow bin.
_MAX_EDITS = 50  # edits-per-read histogram: bins 0.._MAX_EDITS (last = overflow)
_MAX_FRAGLEN = 1000  # fragment-length histogram: bins 0.._MAX_FRAGLEN (overflow)
# Per-read edit-rate histogram: _RATE_BINS equal bins over [0, 1]. Real data piles
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
_SAMPLE_SEED = 20260910


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


class _Stats:
    """Accumulator for QC counts; instances merge with :meth:`update`."""

    def __init__(self) -> None:
        self.total = 0
        self.unmapped = 0
        self.duplicate = 0
        self.secondary = 0
        self.supplementary = 0
        self.proper_pair = 0
        self.passing = 0
        self.total_opportunities = 0
        self.total_edits = 0
        self.edits_per_read = np.zeros(_MAX_EDITS + 1, dtype=np.int64)
        self.fraglen = np.zeros(_MAX_FRAGLEN + 1, dtype=np.int64)
        # Per-read edit rate (edited C/G over editable C/G): distribution + mean.
        self.edit_rate_hist = np.zeros(_RATE_BINS, dtype=np.int64)
        self.edit_rate_sum = 0.0
        self.edit_rate_n = 0
        # Deaminase motif: per-position A/C/G/T counts over the window around
        # each editing event (centre excluded, G-edits reverse-complemented).
        self.motif_pwm = np.zeros((_MOTIF_WINDOW, 4), dtype=np.int64)
        self.motif_events = 0
        # context -> [edits, opportunities]
        self.context: dict[str, list[int]] = defaultdict(lambda: [0, 0])

    def update(self, other: _Stats) -> None:
        self.total += other.total
        self.unmapped += other.unmapped
        self.duplicate += other.duplicate
        self.secondary += other.secondary
        self.supplementary += other.supplementary
        self.proper_pair += other.proper_pair
        self.passing += other.passing
        self.total_opportunities += other.total_opportunities
        self.total_edits += other.total_edits
        self.edits_per_read += other.edits_per_read
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


def _process_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    min_mapq: int,
    min_baseq: int,
    sample_fraction: float = 1.0,
) -> _Stats:
    """Accumulate read, editing, context and fragment-length stats for one chrom.

    When ``sample_fraction`` is below 1, each read is kept with that probability
    and skipped before any analysis, which is where the cost is. The draw is
    seeded from the chromosome name, so a rerun samples the same reads.
    """
    stats = _Stats()
    # A string seed is hashed deterministically by `random`, unlike the builtin
    # hash(), which is salted per process and would break reproducibility.
    rng = random.Random(f"{_SAMPLE_SEED}:{chrom}")
    subsampling = sample_fraction < 1.0

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        ref_seq = fasta.fetch(chrom).upper()
        ref_len = len(ref_seq)

        for read in bam.fetch(chrom):
            if subsampling and rng.random() >= sample_fraction:
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

            if not _passes_filters(read, min_mapq):
                continue
            stats.passing += 1

            # Fragment length from properly-paired read1 only (avoid double count).
            if read.is_proper_pair:
                stats.proper_pair += 1
                if read.is_read1 and read.template_length:
                    flen = min(abs(read.template_length), _MAX_FRAGLEN)
                    stats.fraglen[flen] += 1

            seq = read.query_sequence
            if seq is None:
                continue
            quals = read.query_qualities

            read_edits = 0
            read_editable = 0  # reference C/G covered by this read (strand-agnostic)
            read_edited = 0  # of those, how many show C->T or G->A
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if quals is not None and quals[qpos] < min_baseq:
                    continue
                ref_base = ref_seq[rpos]
                read_base = seq[qpos]

                # Per-read edit rate: every reference C or G is an editable base;
                # a C->T or G->A mismatch is an edit. Counted strand-agnostically
                # and without the flank requirement used for context below.
                if ref_base == "C":
                    read_editable += 1
                    if read_base == "T":
                        read_edited += 1
                elif ref_base == "G":
                    read_editable += 1
                    if read_base == "A":
                        read_edited += 1

                if rpos == 0 or rpos >= ref_len - 1:
                    continue  # need flanking bases for the trinucleotide context

                # Strand-agnostic edit calling (matching bam2bw): a reference C
                # may be edited C->T and a reference G may be edited G->A,
                # regardless of read orientation. The G-centred context is
                # reverse-complemented so both are reported as C->T.
                if ref_base == "C":
                    ctx = ref_seq[rpos - 1 : rpos + 2]
                    is_edit = read_base == "T"
                elif ref_base == "G":
                    ctx = _revcomp(ref_seq[rpos - 1 : rpos + 2])
                    is_edit = read_base == "A"
                else:
                    continue

                if "N" in ctx:
                    continue

                stats.total_opportunities += 1
                slot = stats.context[ctx]
                slot[1] += 1
                if is_edit:
                    stats.total_edits += 1
                    slot[0] += 1
                    read_edits += 1

                    # Deaminase motif: accumulate the reference window around the
                    # edited base (centre excluded), unified to the C->T
                    # orientation by reverse-complementing G-centred windows.
                    half = _MOTIF_WINDOW // 2
                    lo = rpos - half
                    hi = lo + _MOTIF_WINDOW
                    if lo >= 0 and hi <= ref_len:
                        window = ref_seq[lo:hi]
                        if ref_base == "G":
                            window = _revcomp(window)
                        for j, b in enumerate(window):
                            if j == half:
                                continue
                            bi = _BASE_IDX.get(b)
                            if bi is not None:
                                stats.motif_pwm[j, bi] += 1
                        stats.motif_events += 1

            stats.edits_per_read[min(read_edits, _MAX_EDITS)] += 1

            if read_editable > 0:
                rate = read_edited / read_editable
                bin_idx = min(int(rate * _RATE_BINS), _RATE_BINS - 1)
                stats.edit_rate_hist[bin_idx] += 1
                stats.edit_rate_sum += rate
                stats.edit_rate_n += 1

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
) -> Iterator[tuple[str, int, bool]]:
    """Yield ``(chrom, centre, is_minus)`` for TSS whose window fits the contig.

    The TSS is the midpoint of the BED record, so both a 1-bp site and a wider
    interval work. Strand is column 6 (BED6); records without it are treated as
    ``+``, except for the common four-column ``chrom start end strand`` TSS
    files, where column 4 is read as the strand instead.
    """
    with open(tss_path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom not in chrom_sizes:
                continue
            center = (int(parts[1]) + int(parts[2])) // 2
            if center - span < 0 or center + span > chrom_sizes[chrom]:
                continue
            strand_col = 3 if len(parts) == 4 else 5
            is_minus = len(parts) > strand_col and parts[strand_col].strip() == "-"
            yield chrom, center, is_minus


def _tss_enrichment(
    bam_path: str,
    tss_path: str,
    chrom_sizes: dict[str, int],
    min_mapq: int,
    flank: int,
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

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, center, is_minus in _tss_sites(tss_path, chrom_sizes, span):
            start = center - span
            per_tss = np.zeros(n_bins, dtype=np.float64)
            # fetch returns every read overlapping the window, which includes
            # every read whose 5' end falls inside it.
            for read in bam.fetch(chrom, start, center + span):
                if not _passes_filters(read, min_mapq):
                    continue
                # Narrow before the arithmetic: reference_end is optional and
                # `- 1` would run before any check placed after it.
                ref_start, ref_end = read.reference_start, read.reference_end
                if ref_start is None or ref_end is None:
                    continue
                site = ref_end - 1 if read.is_reverse else ref_start
                rel = site - start
                if 0 <= rel < 2 * span:
                    per_tss[rel // bin_size] += 1
            if is_minus:
                per_tss = per_tss[::-1]
            counts += per_tss
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
            float("nan"), positions, counts, mean_per_tss, 0.0, n_tss, bin_size
        )

    profile = mean_per_tss / background
    score = float(profile.max())
    logger.info(
        f"  TSS enrichment = {score:.2f} (peak of the normalised profile over "
        f"{n_tss} TSS, {bin_size}-bp bins, +/-{span} bp)"
    )
    return _TssResult(score, positions, counts, profile, background, n_tss, bin_size)


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
                [f"{pos:g}", int(count), f"{mean:.6g}", f"{float(norm):.6g}"]
            )


def _histogram_summary(hist: np.ndarray) -> dict[str, float]:
    """Mean and median of a value-indexed integer histogram."""
    total = int(hist.sum())
    if total == 0:
        return {"n": 0, "mean": 0.0, "median": 0.0}
    values = np.arange(len(hist))
    mean = float((values * hist).sum() / total)
    cumulative = np.cumsum(hist)
    median = float(np.searchsorted(cumulative, (total + 1) / 2.0))
    return {"n": total, "mean": mean, "median": median}


def _build_metrics(
    stats: _Stats,
    tss: _TssResult | None,
    tss_csv_name: str | None,
) -> dict:
    edit_rate = (
        stats.total_edits / stats.total_opportunities
        if stats.total_opportunities
        else 0.0
    )
    per_read = _histogram_summary(stats.edits_per_read)
    fraglen = _histogram_summary(stats.fraglen)

    # Per-read edit rate (edited C/G over editable C/G). Mean is exact; median is
    # taken from the histogram bin centres.
    rate_mean = stats.edit_rate_sum / stats.edit_rate_n if stats.edit_rate_n else 0.0
    if stats.edit_rate_n:
        cum = np.cumsum(stats.edit_rate_hist)
        med_bin = int(np.searchsorted(cum, (stats.edit_rate_n + 1) / 2.0))
        med_bin = min(med_bin, _RATE_BINS - 1)
        rate_median = (med_bin + 0.5) / _RATE_BINS
    else:
        rate_median = 0.0

    context = {
        ctx: {
            "edits": e,
            "opportunities": o,
            "edit_fraction": (e / o if o else 0.0),
        }
        for ctx, (e, o) in sorted(stats.context.items())
    }

    metrics: dict = {
        "reads": {
            "total": stats.total,
            "passing": stats.passing,
            "unmapped": stats.unmapped,
            "duplicate": stats.duplicate,
            "duplicate_rate": (stats.duplicate / stats.total if stats.total else 0.0),
            "secondary": stats.secondary,
            "supplementary": stats.supplementary,
            "proper_pair": stats.proper_pair,
            "proper_pair_rate": (
                stats.proper_pair / stats.passing if stats.passing else 0.0
            ),
        },
        "editing": {
            "total_opportunities": stats.total_opportunities,
            "total_edits": stats.total_edits,
            "global_edit_rate": edit_rate,
            "mean_edits_per_read": per_read["mean"],
            "median_edits_per_read": per_read["median"],
        },
        "edit_rate_per_read": {
            "n_reads": stats.edit_rate_n,
            "mean": rate_mean,
            "median": rate_median,
            "histogram": [int(c) for c in stats.edit_rate_hist],
            "bin_edges": [round(i / _RATE_BINS, 4) for i in range(_RATE_BINS + 1)],
        },
        "context": context,
        "fragment_length": {
            "mean": fraglen["mean"],
            "median": fraglen["median"],
            "n_pairs": fraglen["n"],
        },
        "motif": {
            "window": _MOTIF_WINDOW,
            "n_events": stats.motif_events,
        },
    }
    if tss is not None:
        # The profile itself is not duplicated here: it is a few hundred numbers
        # and its one home is the CSV, named below so the JSON still points at it.
        metrics["tss_enrichment"] = {
            "score": tss.score,
            "n_tss": tss.n_tss,
            "flank": int(len(tss.counts) // 2 * tss.bin_size),
            "bin_size": tss.bin_size,
            "background": tss.background,
            "total_insertions": int(tss.counts.sum()),
            "profile_csv": tss_csv_name,
        }
    return metrics


def _figure_base64(metrics: dict, stats: _Stats) -> str:
    """Render the multi-panel summary figure and return it as base64 PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, ScalarFormatter

    n_panels = 4
    # One panel per row so each plot is large and readable in the report.
    fig, axes = plt.subplots(n_panels, 1, figsize=(9, 3.4 * n_panels))

    # Panel 1: trinucleotide context edit fraction.
    ctx_items = sorted(
        metrics["context"].items(),
        key=lambda kv: kv[1]["edit_fraction"],
        reverse=True,
    )
    labels = [k for k, _ in ctx_items]
    fracs = [v["edit_fraction"] for _, v in ctx_items]
    axes[0].bar(range(len(labels)), fracs, color="#c0392b")
    axes[0].set_xticks(range(len(labels)))
    axes[0].set_xticklabels(labels, rotation=90, fontsize=6)
    axes[0].set_ylabel("edit fraction")
    axes[0].set_title("Trinucleotide context bias")

    # Panel 2: edits per read (raw count).
    hist = stats.edits_per_read
    axes[1].bar(np.arange(len(hist)), hist, color="#2c7fb8")
    axes[1].set_xlabel("edits per read")
    axes[1].set_ylabel("reads")
    axes[1].set_title("Edits per read")

    # Panel 3: per-read edit rate (edited C/G over editable C/G). Most reads sit
    # below 0.1, so the x axis is symlog: log above 0.01, linear below it so that
    # the (populated) zero-rate bin is still representable. `stairs` is used
    # instead of `bar` because a bar's width is in data units and would be
    # distorted by the non-linear scale.
    rate_hist = stats.edit_rate_hist
    edges = np.arange(len(rate_hist) + 1) / len(rate_hist)
    axes[2].stairs(rate_hist, edges, fill=True, color="#e6550d")
    axes[2].set_xscale("symlog", linthresh=0.01, linscale=0.5)
    axes[2].set_xlim(0, 1)
    axes[2].set_xticks([0, 0.01, 0.1, 1])
    axes[2].xaxis.set_major_formatter(ScalarFormatter())
    axes[2].xaxis.set_minor_locator(NullLocator())
    axes[2].set_xlabel("edit rate per read")
    axes[2].set_ylabel("reads")
    axes[2].set_title(
        f"Per-read edit rate (mean {metrics['edit_rate_per_read']['mean']:.3f})"
    )

    # Panel 4: fragment length.
    fl = stats.fraglen
    axes[3].plot(np.arange(len(fl)), fl, color="#31a354")
    axes[3].set_xlabel("fragment length (bp)")
    axes[3].set_ylabel("pairs")
    axes[3].set_title("Fragment length")

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _tss_figure_base64(tss: _TssResult) -> str:
    """Render the aggregate TSS profile with the score marked on it; base64 PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.plot(tss.positions, tss.profile, color="#756bb1", lw=1.4)
    # Labelled through the legend rather than inline: the flanks sit on this
    # line by construction, so any text next to it lands on the curve.
    ax.axhline(1.0, ls="--", lw=0.8, color="#999", label="flank background")
    ax.legend(loc="upper left", fontsize=9, frameon=False)

    # NaN compares false against itself; that is the "background was zero" case,
    # where there is a profile to look at but no score to mark on it.
    if tss.score == tss.score:
        peak = int(np.argmax(tss.profile))
        ax.axhline(tss.score, ls=":", lw=0.9, color="#d95f02")
        ax.plot(tss.positions[peak], tss.score, "o", color="#d95f02", ms=5)
        label = f"TSS enrichment = {tss.score:.2f}"
    else:
        label = "TSS enrichment = n/a (no flank background)"
    ax.text(
        0.98,
        0.94,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12,
        fontweight="bold",
        color="#d95f02",
        bbox={
            "facecolor": "white",
            "edgecolor": "#e2e5ea",
            "boxstyle": "round,pad=0.4",
        },
    )

    ax.set_xlabel("distance from TSS (bp)")
    ax.set_ylabel("normalized insertions")
    ax.set_title(
        f"TSS enrichment ({tss.n_tss:,} TSS, {tss.bin_size}-bp bins, ENCODE style)"
    )
    ax.margins(x=0)
    # Headroom so the score box clears the peak it annotates.
    top = float(np.nanmax(tss.profile))
    if top > 0:
        ax.set_ylim(top=top * 1.18)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


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


def _motif_logo_base64(pwm: np.ndarray) -> str | None:
    """Render the deaminase motif logo from a count PWM; base64 PNG, or None."""
    if pwm.sum() == 0:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import logomaker  # noqa: E402
    import matplotlib.pyplot as plt  # noqa: E402

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
        "total": "Total read records in the BAM, counted before any filtering.",
        "passing": (
            "Reads kept after filtering (primary, non-duplicate, not QC-fail, "
            "MAPQ &ge; min_mapq). All editing, context and fragment metrics are "
            "computed from these reads only."
        ),
        "unmapped": "Reads flagged unmapped (SAM flag 0x4).",
        "duplicate": "Reads flagged as PCR/optical duplicates (flag 0x400).",
        "duplicate_rate": (
            "Duplicates as a fraction of all reads. High values indicate low "
            "library complexity or over-amplification."
        ),
        "secondary": "Secondary alignments (flag 0x100).",
        "supplementary": "Supplementary / chimeric alignments (flag 0x800).",
        "proper_pair": "Passing reads flagged as a proper pair (flag 0x2).",
        "proper_pair_rate": (
            "Proper pairs as a fraction of passing reads. Low values can signal "
            "insert-size or mapping problems."
        ),
    },
    "editing": {
        "total_opportunities": (
            "Editable reference positions covered by passing reads: any "
            "reference C or G (counted regardless of read orientation), with "
            "both flanking bases present and base quality &ge; min_baseq."
        ),
        "total_edits": (
            "Opportunities showing a deamination event &mdash; a C&rarr;T "
            "mismatch at a reference C or a G&rarr;A mismatch at a reference G, "
            "regardless of read orientation (matching bam2bw)."
        ),
        "global_edit_rate": (
            "total_edits / total_opportunities &mdash; the single most important "
            "signal-quality number. A successful deaminase treatment pushes this "
            "well above the sequencing-error background."
        ),
        "mean_edits_per_read": (
            "Average number of edits per read. Deaminase reads carry many edits, "
            "unlike the two Tn5 cut sites of a standard ATAC read."
        ),
        "median_edits_per_read": "Median number of edits per read.",
    },
    "edit_rate_per_read": {
        "n_reads": "Reads with at least one editable C/G base.",
        "mean": (
            "Mean of the per-read edit rate, where a read's rate = (edited C/G) "
            "/ (editable C/G), counted strand-agnostically over every reference "
            "C and G the read covers. Normalising by the number of editable "
            "bases makes reads comparable regardless of length or composition."
        ),
        "median": "Median of the per-read edit-rate distribution.",
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
            "defined by the ENCODE ATAC-seq pipeline. Higher is better; ENCODE "
            "calls &ge;5 acceptable and &ge;7 ideal for human ATAC, and the "
            "ACCESS-ATAC preprint reports ~13&ndash;14 for a good library."
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


def _fmt(value: object) -> str:
    """Format a metric value for display (floats to 4 sig figs, ints with commas)."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if value != value:  # NaN
            return "n/a"
        return f"{value:.4g}"
    return str(value)


def _html_section(title: str, intro: str, section_key: str, data: dict) -> str:
    """Build one HTML section: heading, intro paragraph, and a metric table."""
    rows = []
    docs = _METRIC_DOCS.get(section_key, {})
    for key, value in data.items():
        desc = docs.get(key, "")
        rows.append(
            f"<tr><td class='k'>{key}</td><td class='v'>{_fmt(value)}</td>"
            f"<td class='d'>{desc}</td></tr>"
        )
    return (
        f"<section><h2>{title}</h2><p class='intro'>{intro}</p>"
        "<table><thead><tr><th>Metric</th><th>Value</th>"
        "<th>Meaning</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _render_html(
    metrics: dict,
    img_b64: str | None,
    motif_b64: str | None,
    tss_b64: str | None,
    bam_path: str,
    fasta_path: str,
    out_name: str,
) -> str:
    """Build a self-contained, MultiQC-style HTML QC report."""
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Highlight cards for the headline numbers.
    cards = [
        ("Passing reads", _fmt(metrics["reads"]["passing"])),
        ("Global edit rate", _fmt(metrics["editing"]["global_edit_rate"])),
        ("Mean edits/read", _fmt(metrics["editing"]["mean_edits_per_read"])),
        ("Mean edit rate/read", _fmt(metrics["edit_rate_per_read"]["mean"])),
        ("Duplicate rate", _fmt(metrics["reads"]["duplicate_rate"])),
    ]
    if "tss_enrichment" in metrics:
        cards.append(("TSS enrichment", _fmt(metrics["tss_enrichment"]["score"])))
    cards_html = "".join(
        f"<div class='card'><div class='cval'>{v}</div>"
        f"<div class='clab'>{lab}</div></div>"
        for lab, v in cards
    )

    img_html = (
        f"<section><h2>Summary plots</h2>"
        f"<img alt='QC summary' src='data:image/png;base64,{img_b64}'></section>"
        if img_b64
        else ""
    )

    motif_html = (
        "<section><h2>Deaminase sequence motif</h2>"
        "<p class='intro'>Sequence logo of the reference bases flanking edited "
        "cytosines, built directly from the editing events in the BAM "
        f"({_fmt(metrics['motif']['n_events'])} events, "
        f"{metrics['motif']['window']}-bp window). The edited base itself is "
        "excluded and G&rarr;A events are reverse-complemented into the "
        "C&rarr;T orientation, so the logo shows the enzyme's flanking-sequence "
        "preference (e.g. DddA's <code>TC</code> bias).</p>"
        f"<img alt='Deaminase motif' src='data:image/png;base64,{motif_b64}'></section>"
        if motif_b64
        else ""
    )

    tss_html = ""
    if tss_b64 and "tss_enrichment" in metrics:
        tss = metrics["tss_enrichment"]
        csv_note = (
            f" The plotted numbers are written to <code>{tss['profile_csv']}</code>."
            if tss.get("profile_csv")
            else ""
        )
        tss_html = (
            "<section><h2>TSS enrichment</h2>"
            "<p class='intro'>Insertion 5' ends aggregated over "
            f"{_fmt(tss['n_tss'])} transcription start sites, following the "
            "ENCODE ATAC-seq pipeline: a &plusmn;"
            f"{_fmt(tss['flank'])}&nbsp;bp window binned at "
            f"{_fmt(tss['bin_size'])}&nbsp;bp, minus-strand TSS flipped so "
            "upstream is always on the left, then divided by the mean signal in "
            "the outermost 100&nbsp;bp on each side so the flanks sit at 1. The "
            "score is the peak of that normalised profile."
            f"{csv_note}</p>"
            f"<img alt='TSS enrichment' src='data:image/png;base64,{tss_b64}'>"
            "</section>"
        )

    sections = [
        _html_section(
            "Read statistics",
            "Composition of the alignment file before and after filtering.",
            "reads",
            metrics["reads"],
        ),
        _html_section(
            "Editing statistics",
            "Strand-agnostic deamination signal: how many editable C/G "
            "positions were seen and how many were edited (C&rarr;T or "
            "G&rarr;A, matching bam2bw).",
            "editing",
            metrics["editing"],
        ),
        _html_section(
            "Per-read edit rate",
            "Distribution of each read's edited-fraction of its editable C/G "
            "bases. Plotted as its own panel in the summary figure above.",
            "edit_rate_per_read",
            {
                k: v
                for k, v in metrics["edit_rate_per_read"].items()
                if k not in ("histogram", "bin_edges")
            },
        ),
        _html_section(
            "Fragment length",
            "Insert-size distribution from properly-paired reads.",
            "fragment_length",
            metrics["fragment_length"],
        ),
    ]
    if "tss_enrichment" in metrics:
        sections.append(
            _html_section(
                "TSS enrichment metrics",
                "The numbers behind the profile plotted above.",
                "tss_enrichment",
                metrics["tss_enrichment"],
            )
        )

    # Trinucleotide context gets a bespoke table (one row per context).
    ctx_rows = "".join(
        f"<tr><td class='k'>{ctx}</td><td class='v'>{_fmt(d['edit_fraction'])}</td>"
        f"<td class='v'>{_fmt(d['edits'])}</td>"
        f"<td class='v'>{_fmt(d['opportunities'])}</td></tr>"
        for ctx, d in sorted(
            metrics["context"].items(),
            key=lambda kv: kv[1]["edit_fraction"],
            reverse=True,
        )
    )
    ctx_section = (
        "<section><h2>Trinucleotide context bias</h2>"
        "<p class='intro'>Edit fraction per cytosine-centred trinucleotide "
        "(G&rarr;A events are reverse-complemented into the "
        "C&rarr;T orientation, so both strands are unified). This is the "
        "enzyme's sequence-preference fingerprint: DddA strongly prefers "
        "<code>TC</code> contexts, while relaxed-bias enzymes (DddSs / "
        "SsdAtox) edit more uniformly. A strongly skewed profile means "
        "downstream footprinting should apply enzyme-bias correction.</p>"
        "<table><thead><tr><th>Context</th><th>Edit fraction</th>"
        "<th>Edits</th><th>Opportunities</th></tr></thead><tbody>"
        + ctx_rows
        + "</tbody></table></section>"
    )

    style = """
    body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
      color:#222;margin:0;background:#f5f6f8;}
    .container{max-width:1000px;margin:0 auto;padding:24px;}
    h1{color:#16767a;margin-bottom:2px;}
    .meta{color:#666;font-size:13px;margin-top:0;}
    .cards{display:flex;flex-wrap:wrap;gap:12px;margin:20px 0;}
    .card{background:#fff;border:1px solid #e2e5ea;border-radius:8px;
      padding:14px 18px;min-width:140px;flex:1;text-align:center;
      box-shadow:0 1px 2px rgba(0,0,0,.04);}
    .cval{font-size:22px;font-weight:700;color:#16767a;}
    .clab{font-size:12px;color:#666;margin-top:4px;}
    section{background:#fff;border:1px solid #e2e5ea;border-radius:8px;
      padding:8px 20px 18px;margin:18px 0;}
    h2{color:#16767a;font-size:18px;border-bottom:1px solid #eee;
      padding-bottom:6px;}
    .intro{color:#555;font-size:14px;}
    table{border-collapse:collapse;width:100%;font-size:13px;}
    th,td{text-align:left;padding:6px 10px;border-bottom:1px solid #eee;
      vertical-align:top;}
    th{background:#fafbfc;color:#444;}
    td.k{font-family:monospace;white-space:nowrap;color:#16767a;}
    td.v{font-family:monospace;white-space:nowrap;}
    td.d{color:#555;}
    img{max-width:100%;height:auto;}
    code{background:#f0f2f4;padding:1px 4px;border-radius:3px;}
    footer{color:#999;font-size:12px;text-align:center;margin:24px 0;}
    """

    return (
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>DeamTools QC Report</title>"
        f"<style>{style}</style></head><body><div class='container'>"
        "<h1>DeamTools QC Report</h1>"
        f"<p class='meta'>Sample: <b>{out_name}</b> &middot; Generated "
        f"{generated}<br>BAM: <code>{bam_path}</code><br>"
        f"FASTA: <code>{fasta_path}</code></p>"
        f"<div class='cards'>{cards_html}</div>"
        f"{img_html}"
        f"{tss_html}"
        f"{motif_html}"
        + "".join(sections)
        + ctx_section
        + "<footer>Generated by deamtools qc</footer>"
        "</div></body></html>"
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
) -> dict:
    """Compute QC metrics for a deaminase chromatin-accessibility BAM.

    Writes a machine-readable ``<out_dir>/<out_name>.json`` and a
    self-contained, MultiQC-style ``<out_dir>/<out_name>.html`` report that
    embeds the summary figure and documents every metric inline. With
    ``tss_path``, the aggregate TSS profile behind the report's plot is also
    written to ``<out_dir>/<out_name>.tss_enrichment.csv``.

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
        Number of worker threads for per-chromosome processing.
    tss_flank : int, default 2000
        Half-width (bp) of the window around each TSS for enrichment. The
        ENCODE default is 2000; it is rounded down to a whole number of 10-bp
        bins.
    plot : bool, default True
        Whether to render and embed the summary figure in the HTML report.
    n_reads : int, optional
        Subsample to approximately this many reads instead of using all of
        them, for a faster pass over a large BAM. Reads are drawn uniformly
        across the genome (each kept with probability ``n_reads / total``,
        which is read from the BAM index, so nothing is scanned twice), and
        the draw is seeded, so a rerun samples the same reads. Ignored when
        the BAM already holds fewer reads than this.

        Rates and distributions -- editing rate, duplicate rate, context bias,
        the motif PWM, fragment lengths -- are unbiased under this sampling.
        The absolute counts in the report are counts *of the sample*, not
        estimates of the whole file, and the ``sampling`` block of the metrics
        records the fraction so they can be scaled if needed. TSS enrichment
        always uses every read in its windows, since it is a ratio computed
        over a small part of the genome and subsampling it would only add
        noise.

    Returns
    -------
    dict
        The metrics dictionary (also written to ``<out_dir>/<out_name>.json``).
    """
    logger.info("Running qc")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_sizes = get_chrom_sizes_from_bam(bam)
        # The index carries per-contig counts, so the total costs no scan.
        total_reads = sum(st.mapped + st.unmapped for st in bam.get_index_statistics())
    chroms = list(chrom_sizes.keys())

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

    logger.info(f"Processing {len(chroms)} chromosome(s) with {threads} thread(s)")

    merged = _Stats()
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(
                _process_chrom,
                bam_path,
                fasta_path,
                c,
                min_mapq,
                min_baseq,
                sample_fraction,
            ): c
            for c in chroms
        }
        for future in as_completed(futures):
            merged.update(future.result())

    logger.info(
        f"  {merged.passing} passing read(s); "
        f"{merged.total_edits}/{merged.total_opportunities} edits/opportunities"
    )

    tss: _TssResult | None = None
    if tss_path is not None:
        logger.info(f"TSS:   {tss_path}")
        tss = _tss_enrichment(bam_path, tss_path, chrom_sizes, min_mapq, tss_flank)

    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f"{out_name}.json")
    html_path = os.path.join(out_dir, f"{out_name}.html")
    tss_csv_name = f"{out_name}.tss_enrichment.csv" if tss is not None else None

    metrics = _build_metrics(merged, tss, tss_csv_name)
    metrics["sampling"] = {
        "subsampled": sample_fraction < 1.0,
        "requested_reads": n_reads,
        "fraction": sample_fraction,
        "reads_in_bam": total_reads,
    }

    if tss is not None and tss_csv_name is not None:
        tss_csv_path = os.path.join(out_dir, tss_csv_name)
        logger.info(f"Writing {tss_csv_path}")
        _write_tss_csv(tss_csv_path, tss)

    logger.info(f"Writing {json_path}")
    with open(json_path, "w") as f:
        json.dump(metrics, f, indent=2)

    img_b64 = _figure_base64(metrics, merged) if plot else None
    motif_b64 = _motif_logo_base64(merged.motif_pwm) if plot else None
    tss_b64 = _tss_figure_base64(tss) if plot and tss is not None else None

    logger.info(f"Writing {html_path}")
    with open(html_path, "w") as f:
        f.write(
            _render_html(
                metrics, img_b64, motif_b64, tss_b64, bam_path, fasta_path, out_name
            )
        )

    logger.info("Done")
    return metrics
