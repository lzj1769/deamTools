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
* **TSS enrichment** *(optional)* — the classic ATAC-style enrichment of Tn5
  insertion sites around transcription start sites, computed when a TSS BED is
  supplied.

Results are written as machine-readable JSON plus a self-contained,
MultiQC-style HTML report (``<out_dir>/<out_name>.json`` and ``.html``). The
HTML embeds the multi-panel summary figure and documents the meaning of every
metric inline.
"""

from __future__ import annotations

import base64
import io
import json
import logging
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import numpy as np
import pysam

from deamtools.utils import get_chrom_sizes_from_bam

logger = logging.getLogger(__name__)

# Histograms are stored as fixed-length arrays with a final overflow bin.
_MAX_EDITS = 50  # edits-per-read histogram: bins 0.._MAX_EDITS (last = overflow)
_MAX_FRAGLEN = 1000  # fragment-length histogram: bins 0.._MAX_FRAGLEN (overflow)
_RATE_BINS = 50  # per-read edit-rate histogram: _RATE_BINS equal bins over [0, 1]

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


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
) -> _Stats:
    """Accumulate read, editing, context and fragment-length stats for one chrom."""
    stats = _Stats()
    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        ref_seq = fasta.fetch(chrom).upper()
        ref_len = len(ref_seq)

        for read in bam.fetch(chrom):
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

            stats.edits_per_read[min(read_edits, _MAX_EDITS)] += 1

            if read_editable > 0:
                rate = read_edited / read_editable
                bin_idx = min(int(rate * _RATE_BINS), _RATE_BINS - 1)
                stats.edit_rate_hist[bin_idx] += 1
                stats.edit_rate_sum += rate
                stats.edit_rate_n += 1

    return stats


def _tss_enrichment(
    bam_path: str,
    tss_path: str,
    chrom_sizes: dict[str, int],
    min_mapq: int,
    flank: int,
) -> tuple[float, np.ndarray]:
    """ATAC-style TSS enrichment from Tn5 insertion sites.

    The insertion site of a read is its 5' end (``reference_start`` for forward
    reads, ``reference_end - 1`` for reverse reads). Insertions are aggregated
    into a ``2 * flank + 1`` profile centred on every TSS, normalised by the
    mean insertion density in the outer flanks, and the enrichment score is the
    mean of the normalised profile in the central ``+/-50`` bp window.

    Returns the enrichment score and the normalised profile.
    """
    width = 2 * flank + 1
    profile = np.zeros(width, dtype=np.float64)
    n_tss = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with open(tss_path) as f:
            for line in f:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split("\t")
                chrom = parts[0]
                if chrom not in chrom_sizes:
                    continue
                center = (int(parts[1]) + int(parts[2])) // 2
                start = center - flank
                end = center + flank + 1
                if start < 0 or end > chrom_sizes[chrom]:
                    continue
                n_tss += 1
                for read in bam.fetch(chrom, start, end):
                    if not _passes_filters(read, min_mapq):
                        continue
                    site = (
                        read.reference_end - 1
                        if read.is_reverse
                        else read.reference_start
                    )
                    if site is None:
                        continue
                    rel = site - start
                    if 0 <= rel < width:
                        profile[rel] += 1

    if n_tss == 0:
        logger.warning("  no usable TSS found; skipping TSS enrichment")
        return float("nan"), profile

    # Background = mean insertion density in the outermost 100 bp on each side.
    edge = min(100, flank)
    background = np.concatenate([profile[:edge], profile[-edge:]]).mean()
    if background <= 0:
        logger.warning("  zero TSS flank background; enrichment undefined")
        return float("nan"), profile
    normalized = profile / background
    half = min(50, flank)
    score = float(normalized[flank - half : flank + half + 1].mean())
    logger.info(f"  TSS enrichment = {score:.2f} (over {n_tss} TSS)")
    return score, normalized


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
    tss_score: float | None,
    tss_profile: np.ndarray | None,
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
    rate_mean = (
        stats.edit_rate_sum / stats.edit_rate_n if stats.edit_rate_n else 0.0
    )
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
    }
    if tss_score is not None:
        metrics["tss_enrichment"] = {
            "score": tss_score,
            "profile": [round(float(v), 4) for v in tss_profile],
        }
    return metrics


def _figure_base64(metrics: dict, stats: _Stats) -> str:
    """Render the multi-panel summary figure and return it as base64 PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_tss = "tss_enrichment" in metrics
    n_panels = 5 if has_tss else 4
    # One panel per row so each plot is large and readable in the report.
    fig, axes = plt.subplots(
        n_panels, 1, figsize=(9, 3.4 * n_panels)
    )

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

    # Panel 3: per-read edit rate (edited C/G over editable C/G).
    rate_hist = stats.edit_rate_hist
    n_bins = len(rate_hist)
    centers = (np.arange(n_bins) + 0.5) / n_bins
    axes[2].bar(centers, rate_hist, width=1.0 / n_bins, color="#e6550d")
    axes[2].set_xlim(0, 1)
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

    # Panel 5: TSS enrichment.
    if has_tss:
        prof = metrics["tss_enrichment"]["profile"]
        flank = (len(prof) - 1) // 2
        x = np.arange(-flank, flank + 1)
        axes[4].plot(x, prof, color="#756bb1")
        axes[4].axhline(1.0, ls="--", lw=0.8, color="grey")
        axes[4].set_xlabel("distance from TSS (bp)")
        axes[4].set_ylabel("normalized insertions")
        axes[4].set_title(
            f"TSS enrichment = {metrics['tss_enrichment']['score']:.2f}"
        )

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
            "ATAC-style enrichment of Tn5 insertion 5' ends in the central "
            "&plusmn;50 bp around TSS, normalised to the flank background. Higher "
            "is better; &gt;6&ndash;10 indicates good accessible-chromatin "
            "enrichment."
        ),
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
                "TSS enrichment",
                "Accessible-chromatin enrichment of Tn5 insertions around "
                "transcription start sites.",
                "tss_enrichment",
                {"score": metrics["tss_enrichment"]["score"]},
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
) -> dict:
    """Compute QC metrics for a deaminase chromatin-accessibility BAM.

    Writes two files: a machine-readable ``<out_dir>/<out_name>.json`` and a
    self-contained, MultiQC-style ``<out_dir>/<out_name>.html`` report that
    embeds the summary figure and documents every metric inline.

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
        BED file of transcription start sites. When supplied, an ATAC-style TSS
        enrichment score and profile are computed.
    min_mapq : int, default 20
        Minimum read mapping quality.
    min_baseq : int, default 20
        Minimum base quality for a position to count as an editing opportunity.
    threads : int, default 1
        Number of worker threads for per-chromosome processing.
    tss_flank : int, default 2000
        Half-width (bp) of the window around each TSS for enrichment.
    plot : bool, default True
        Whether to render and embed the summary figure in the HTML report.

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
    chroms = list(chrom_sizes.keys())
    logger.info(f"Processing {len(chroms)} chromosome(s) with {threads} thread(s)")

    merged = _Stats()
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(
                _process_chrom, bam_path, fasta_path, c, min_mapq, min_baseq
            ): c
            for c in chroms
        }
        for future in as_completed(futures):
            merged.update(future.result())

    logger.info(
        f"  {merged.passing} passing read(s); "
        f"{merged.total_edits}/{merged.total_opportunities} edits/opportunities"
    )

    tss_score: float | None = None
    tss_profile: np.ndarray | None = None
    if tss_path is not None:
        logger.info(f"TSS:   {tss_path}")
        tss_score, tss_profile = _tss_enrichment(
            bam_path, tss_path, chrom_sizes, min_mapq, tss_flank
        )

    metrics = _build_metrics(merged, tss_score, tss_profile)

    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f"{out_name}.json")
    html_path = os.path.join(out_dir, f"{out_name}.html")

    logger.info(f"Writing {json_path}")
    with open(json_path, "w") as f:
        json.dump(metrics, f, indent=2)

    img_b64 = _figure_base64(metrics, merged) if plot else None

    logger.info(f"Writing {html_path}")
    with open(html_path, "w") as f:
        f.write(_render_html(metrics, img_b64, bam_path, fasta_path, out_name))

    logger.info("Done")
    return metrics
