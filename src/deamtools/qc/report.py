"""Offline QC dashboard rendering, independent of BAM processing.

Interpretation is descriptive: no assay-specific quality cutoffs are invented.
The overall status covers evaluability and arithmetic consistency only.
"""

from __future__ import annotations

import base64
import copy
import csv
import html
import json
import math
import re
import shlex
import shutil
from pathlib import Path
from urllib.parse import quote

from deamtools.qc.report_assets import CSS, JS

_PARAMETER_DOCS = {
    "min_mapq": "Minimum MAPQ required for a mapped read to enter downstream QC calculations.",
    "min_baseq": "Minimum base quality required when counting editable C/G opportunities and editing events.",
    "tss_flank": "Half-width in bp around each TSS; rounded down to whole profile bins.",
    "n_reads": "Approximate target number of analyzed read records; deterministic sampling keeps mates together. None means no subsampling.",
    "threads": "Number of worker processes for chromosome-level QC accumulation.",
    "plot": "Whether the report includes plots; metric and CSV outputs are retained either way.",
    "logo_scale": "Frequency or information-content (bits) scale for the sequence logo.",
    "report_dir": "Portable report directory containing report.html, assets and companion data.",
    "redact_paths": "Remove directory components from report provenance and recorded command paths; local QC JSON retains full provenance.",
}


def _path_basename(value: str) -> str:
    """Remove directory components of Unix, Windows and relative path strings."""
    return value.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1] or "[root]"


def _redacted_command(command: str) -> str:
    try:
        tokens = shlex.split(command)
    except ValueError:
        # A malformed command must not fall back to leaking the original text.
        return "[command omitted: could not safely parse recorded arguments]"
    result = []
    for token in tokens:
        if token.startswith("-") and "=" in token:
            flag, value = token.split("=", 1)
            result.append(
                f"{flag}={_path_basename(value)}"
                if "/" in value or "\\" in value
                else token
            )
        else:
            result.append(
                _path_basename(token) if "/" in token or "\\" in token else token
            )
    return shlex.join(result)


def _shareable_metrics(metrics: dict) -> dict:
    """Redact only a rendering copy; never mutate measurements or local provenance."""
    result = copy.deepcopy(metrics)

    def sanitize(value, field=""):
        if isinstance(value, dict):
            return {
                key: (
                    _redacted_command(item)
                    if key == "command" and isinstance(item, str)
                    else sanitize(item, key)
                )
                for key, item in value.items()
            }
        if isinstance(value, list):
            return [sanitize(item, field) for item in value]
        if isinstance(value, str) and (
            field in {"bam", "fasta", "tss", "report_dir", "out_dir"}
            or field.endswith(("_path", "_dir", "_file"))
            or value.startswith(("/", "~/", "../", "./", "file://", "\\"))
            or re.match(r"^[A-Za-z]:[\\/]", value)
        ):
            return _path_basename(value)
        return value

    if "provenance" in result:
        result["provenance"] = sanitize(result["provenance"])
    return result


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def _number(value: object, *, percent: bool = False, compact: bool = False) -> str:
    if not _finite(value):
        return "N/A"
    number = float(value)  # type: ignore[arg-type]
    if percent:
        return f"{number * 100:.2f}%"
    if compact:
        for divisor, suffix in ((1e9, "B"), (1e6, "M"), (1e3, "K")):
            if abs(number) >= divisor:
                return f"{number / divisor:.2f} {suffix}"
    return f"{int(number):,}" if number.is_integer() else f"{number:,.2f}"


def _plain(value: object) -> str:
    if value is None or isinstance(value, float) and not math.isfinite(value):
        return "N/A"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, (int, float)):
        return _number(value)
    if isinstance(value, dict):
        return (
            "; ".join(f"{key}: {_plain(val)}" for key, val in value.items())
            or "None recorded"
        )
    return str(value)


def _info(description: str) -> str:
    # Metric descriptions may contain legacy HTML entities; tooltips are plain text.
    text = html.unescape(re.sub(r"<[^>]*>", "", description))
    return (
        f"<button type='button' class='info' title='{_escape(text)}' aria-label='{_escape(text)}'>ⓘ"
        f"<span class='tip' role='tooltip'>{_escape(text)}</span></button>"
    )


def _table(
    data: dict, descriptions: dict | None = None, *, rate_keys: tuple = ()
) -> str:
    rows = []
    for key, value in data.items():
        rate = key in rate_keys or key.endswith(("_rate", "_fraction"))
        shown = (
            _number(value, percent=True) if rate and _finite(value) else _plain(value)
        )
        explanation = (descriptions or {}).get(key, "")
        rows.append(
            f"<tr><td>{_escape(key)}</td><td class='value' title='{_escape(value)}'>{_escape(shown)}</td><td>{_info(explanation) if explanation else ""}</td></tr>"
        )
    return (
        "<div class='table-scroll'><table><thead><tr><th>Metric</th><th>Value</th><th>Info</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
    )


def _details(title: str, body: str, *, identity: str = "") -> str:
    anchor = f" id='{identity}'" if identity else ""
    return f"<details{anchor}><summary>{_escape(title)}</summary>{body}</details>"


def _bias_summary(metrics: dict) -> dict:
    contexts = {
        key: value
        for key, value in metrics.get("context", {}).items()
        if value.get("opportunities", 0) > 0 and _finite(value.get("edit_fraction"))
    }
    if not contexts:
        return {
            "top": None,
            "rate": None,
            "baseline": None,
            "fold": None,
            "varies": False,
        }
    top = max(contexts, key=lambda key: contexts[key]["edit_fraction"])
    opportunities = sum(value["opportunities"] for value in contexts.values())
    baseline = sum(value["edits"] for value in contexts.values()) / opportunities
    rate = contexts[top]["edit_fraction"]
    return {
        "top": top,
        "rate": rate,
        "baseline": baseline,
        "fold": rate / baseline if baseline else None,
        "varies": len({value["edit_fraction"] for value in contexts.values()}) > 1,
    }


def assess_qc(metrics: dict) -> dict:
    """Assess data availability/consistency and generate descriptive advisories.

    PASS is not an assay-quality certification. Checks use logical invariants,
    never empirical edit-rate, duplication or TSS thresholds.
    """
    reads, editing = metrics.get("reads", {}), metrics.get("editing", {})
    total, passing = reads.get("total"), reads.get("passing")
    opportunities, edits = editing.get("total_opportunities"), editing.get(
        "total_edits"
    )
    checks = [
        ("Passing reads available", _finite(passing) and passing > 0),
        ("Editable bases available", _finite(opportunities) and opportunities > 0),
        (
            "Read counts consistent",
            _finite(total) and _finite(passing) and 0 <= passing <= total,
        ),
        (
            "Edit counts consistent",
            _finite(opportunities) and _finite(edits) and 0 <= edits <= opportunities,
        ),
    ]
    warnings = []
    recommendations = []
    rate = editing.get("global_edit_rate")
    summary = []
    if _finite(rate):
        summary.append("Editing signal was successfully quantified.")
        if edits == 0:
            warnings.append(
                "No C→T/G→A editing events were observed among the evaluated opportunities."
            )
        recommendations.append(
            "Compare editing activity against matched enzyme, treatment and library controls; no universal edit-rate cutoff is applied."
        )
    else:
        summary.append("Editing activity cannot be evaluated from the available data.")
    duplicate_rate = reads.get("duplicate_rate")
    if duplicate_rate == 0:
        recommendations.append(
            "No duplicate-flagged records were present in the analyzed BAM. This may reflect upstream deduplication."
        )
    elif _finite(duplicate_rate):
        recommendations.append(
            f"Duplicate-flagged records represent {_number(duplicate_rate, percent=True)} of the analyzed records; interpret this alongside library preparation and sequencing depth."
        )
    tss = metrics.get("tss_enrichment", {})
    score = tss.get("score")
    tss_status = metrics.get("status", {}).get("tss")
    if _finite(score):
        summary.append("TSS enrichment was successfully quantified.")
        recommendations.append(
            "Compare TSS enrichment only with the same annotation, window, binning and cut-site definition. This report does not assign ATAC pass/fail thresholds to deaminase libraries."
        )
    elif tss_status not in (None, "not_requested") or tss:
        warnings.append(
            f"TSS enrichment is unavailable ({tss_status or tss.get('status', 'undefined')}); inspect annotation coverage and flank counts."
        )
    annotation = metrics.get("status", {}).get("tss_annotation", {})
    excluded = annotation.get("unknown_contig", 0) + annotation.get("outside_contig", 0)
    if excluded:
        warnings.append(
            f"{excluded:,} TSS annotation records were excluded because of contig names or window boundaries."
        )
    bias = _bias_summary(metrics)
    if bias["varies"]:
        warnings.append(
            "Observed edit rates vary across sequence contexts; this is descriptive, not a statistical test of enzyme bias."
        )
        recommendations.append(
            "Consider context-bias correction before footprint interpretation; distinguish sequence preference from chromatin protection."
        )
    if metrics.get("sampling", {}).get("subsampled"):
        recommendations.append(
            "Read and editing counts describe the sampled records; file_reads is unsampled and TSS uses all reads in its windows."
        )
    recommendations.append(
        "Zero-edit fragments alone do not establish closed chromatin; inspect their editable-base opportunity counts."
    )
    if _finite(rate) and _finite(score):
        summary = ["Editing signal and TSS enrichment were successfully quantified."]
    failed = [name for name, passed in checks if not passed]
    overall = "FAIL" if failed else "PASS"
    return {
        "overall": overall,
        "checks": checks,
        "passed": sum(passed for _, passed in checks),
        "warnings": warnings,
        "failed": failed,
        "summary": " ".join(summary)
        + " See the metric cards and detailed sections below.",
        "recommendations": recommendations,
        "bias": bias,
    }


class _Assets:
    def __init__(self, source: Path, destination: Path | None):
        self.source = source
        self.destination = destination
        self.files: dict[str, str] = {}
        if destination:
            (destination / "assets").mkdir(parents=True, exist_ok=True)

    def csv(self, filename: str | None) -> tuple[list[dict], str]:
        if not filename:
            return [], ""
        # Reports refer only to companion files in their output directory.
        if Path(filename).name != filename:
            return [], ""
        path = self.source / filename
        if not path.is_file():
            return [], ""
        payload = path.read_bytes()
        if self.destination:
            target = self.destination / filename
            if target.resolve() != path.resolve():
                shutil.copyfile(path, target)
            href = quote(filename)
        else:
            href = "data:text/csv;base64," + base64.b64encode(payload).decode("ascii")
        self.files[filename] = href
        with path.open(newline="") as stream:
            return list(csv.DictReader(stream)), href

    def image(self, encoded: str, name: str) -> str:
        if self.destination:
            (self.destination / "assets" / f"{name}.png").write_bytes(
                base64.b64decode(encoded)
            )
            return f"assets/{name}.png"
        return "data:image/png;base64," + encoded

    def svg(self, svg: str, name: str) -> None:
        if self.destination:
            (self.destination / "assets" / f"{name}.svg").write_text(
                svg, encoding="utf-8"
            )

    def download(self, filename: str | None) -> str:
        if filename not in self.files:
            return ""
        return f"<a class='download' href='{self.files[filename]}' download='{_escape(filename)}'>↓ CSV</a>"


def _chart(
    identity: str,
    title: str,
    xs: list[float],
    ys: list[float | None],
    assets: _Assets,
    *,
    xlabel: str,
    ylabel: str,
    labels: list[str] | None = None,
    percent_x: bool = False,
    line: bool = False,
    details: list[str] | None = None,
    csv_name: str | None = None,
) -> str:
    """Render an independent SVG panel with exact hover values and offline controls."""
    valid = [
        (x, float(y), i)
        for i, (x, y) in enumerate(zip(xs, ys, strict=True))
        if y is not None and math.isfinite(y)
    ]
    if not valid:
        return f"<div class='chart'><h3>{_escape(title)}</h3><p class='empty'>No evaluable data for this panel.</p>{assets.download(csv_name)}</div>"
    lo, hi = min(xs), max(xs)
    if hi <= lo:
        hi = lo + 1

    def coordinate(value):
        return (
            value
            if not percent_x or value <= 0.01
            else 0.01 * (1 + math.log(value / 0.01))
        )

    def inverse(value):
        return (
            value
            if not percent_x or value <= 0.01
            else 0.01 * math.exp(value / 0.01 - 1)
        )

    lo, hi = coordinate(lo), coordinate(hi)
    # Category positions sit inside the bounds rather than on the frame edge.
    if labels:
        lo, hi = -0.5, len(labels) - 0.5
    peak = max(float(y) for _, y, _ in valid)
    ymax = peak * 1.08 if peak > 0 else 1.0
    shapes = []
    points = []
    bar_width = min(28, 540 / max(len(xs), 1) * 0.85)
    for x, y, i in valid:
        px, py = (coordinate(x) - lo) / (hi - lo) * 540, 220 - float(y) / ymax * 220
        description = (
            details[i]
            if details
            else f"{xlabel}: {labels[i] if labels else x}; {ylabel}: {y}"
        )
        tip = _escape(description)
        if line:
            points.append(f"{px:.3f},{py:.3f}")
            shapes.append(
                f"<circle class='plot-point' cx='{px:.3f}' cy='{py:.3f}' r='3' fill='#16767a' tabindex='0' aria-label='{tip}'><title>{tip}</title></circle>"
            )
        else:
            shapes.append(
                f"<rect class='plot-point' x='{px-bar_width/2:.3f}' y='{py:.3f}' width='{bar_width:.3f}' height='{max(.5,220-py):.3f}' fill='#16767a' tabindex='0' aria-label='{tip}'><title>{tip}</title></rect>"
            )
    if line:
        shapes.insert(
            0,
            f"<polyline points='{' '.join(points)}' fill='none' stroke='#16767a' stroke-width='1.5' vector-effect='non-scaling-stroke'/>",
        )
    axis = []
    for i in range(5):
        y = 245 - i * 55
        axis.append(
            f"<line x1='65' x2='605' y1='{y}' y2='{y}' stroke='#e3ebeb'/><text x='57' y='{y+4}' text-anchor='end' font-size='11'>{_number(ymax*i/4, compact=True)}</text>"
        )
        val = lo + (hi - lo) * i / 4
        tick = (
            labels[min(len(labels) - 1, max(0, round(val)))]
            if labels
            else _number(inverse(val), percent=percent_x)
        )
        axis.append(
            f"<text data-tick='{i}' x='{65+135*i}' y='266' text-anchor='middle' font-size='11'>{_escape(tick)}</text>"
        )
    svg = (
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 640 310' role='img' aria-label='{_escape(title)}' font-family='system-ui,sans-serif' fill='#36565a'>"
        f"<title>{_escape(title)}</title>{''.join(axis)}"
        f"<svg x='65' y='25' width='540' height='220' viewBox='0 0 540 220' preserveAspectRatio='none' data-viewport='true' overflow='hidden'>{''.join(shapes)}</svg>"
        f"<text x='335' y='298' text-anchor='middle' font-size='12'>{_escape(xlabel)}</text>"
        f"<text x='12' y='16' font-size='11'>{_escape(ylabel)}</text></svg>"
    )
    assets.svg(svg, identity)
    return (
        f"<div class='chart' data-chart='{identity}' data-min='{lo}' data-max='{hi}' data-scale='{'symlog-rate' if percent_x else 'number'}' data-labels='{_escape(json.dumps(labels or []))}'>"
        f"<h3>{_escape(title)}</h3>{svg}<div class='chart-tools'>"
        f"<button type='button' data-zoom='in' aria-label='Zoom in: {_escape(title)}'>＋</button>"
        f"<button type='button' data-zoom='out' aria-label='Zoom out: {_escape(title)}'>−</button>"
        "<button type='button' data-zoom='reset'>Reset</button><span data-zoom-label>1×</span>"
        f"<label>Pan <input aria-label='Pan: {_escape(title)}' type='range' min='0' max='100' value='0'></label>"
        f"<button type='button' data-export>↓ SVG</button>{assets.download(csv_name)}</div></div>"
    )


def write_report(
    metrics: dict,
    source_dir: str,
    out_name: str,
    *,
    plot: bool = True,
    motif_b64: str | None = None,
    report_dir: str | None = None,
    descriptions: dict | None = None,
    html_path: str | None = None,
    redact_paths: bool = False,
) -> str:
    """Write an offline QC dashboard from metrics and companion CSV files.

    Parameters
    ----------
    metrics : dict
        Actual QC measurements and optional recorded provenance.
    source_dir : str
        Directory containing the CSV files referenced by metrics.
    out_name : str
        Sample name fallback and default HTML basename.
    plot : bool
        Include independent SVG panels and a motif logo when available.
    motif_b64 : str, optional
        Existing PNG logo. If omitted, rebuild it from the recorded motif PFM.
    report_dir : str, optional
        Write report.html, external image/vector assets and companion CSVs here.
    descriptions : dict, optional
        Metric definitions for concise tables and accessible tooltips.
    html_path : str, optional
        Override the standalone HTML destination (not used in directory mode).

    redact_paths : bool
        Redact provenance directory components and command paths from the HTML
        and the portable metrics.json copy; source metrics remain unchanged.

    Returns
    -------
    str
        Path to the generated report.
    """
    if redact_paths:
        metrics = _shareable_metrics(metrics)
    source = Path(source_dir)
    destination = Path(report_dir) if report_dir else None
    assets = _Assets(source, destination)
    docs = descriptions or {}
    provenance = metrics.get("provenance", {})
    sample = provenance.get("sample", out_name)
    layout = metrics.get("library_layout", "Not recorded")
    reads, editing = metrics.get("reads", {}), metrics.get("editing", {})
    rates, motif = metrics.get("edit_rate_per_fragment", {}), metrics.get("motif", {})
    tss = metrics.get("tss_enrichment", {})
    assessment = assess_qc(metrics)
    bias = assessment["bias"]
    csvs: dict[str, list[dict]] = {}
    for section in (
        editing,
        rates,
        motif,
        tss,
        metrics.get("fragment_length", {}),
        metrics.get("context_summary", {}),
    ):
        for key, value in section.items():
            if (key.endswith("_csv") or key == "csv") and isinstance(value, str):
                csvs[value], _ = assets.csv(value)
    if destination:
        (destination / "metrics.json").write_text(
            json.dumps(metrics, indent=2, allow_nan=False), encoding="utf-8"
        )

    def table(name, data=None, rate_keys=()):
        return _table(
            metrics.get(name, {}) if data is None else data,
            docs.get(name),
            rate_keys=rate_keys,
        )

    def panel(
        key,
        title,
        section,
        csv_key,
        xkey,
        ykey,
        *,
        xlabel,
        ylabel,
        line=False,
        percent_x=False,
    ):
        filename = section.get(csv_key)
        rows = csvs.get(filename, [])
        if not plot:
            return ""
        xs = [float(row[xkey]) for row in rows]
        ys = [float(row[ykey]) if row.get(ykey) else None for row in rows]
        if key == "edits-distribution" and rows:
            occupied = [i for i, val in enumerate(ys) if val]
            stop = occupied[-1] + 1 if occupied else 1
            xs, ys, rows = xs[:stop], ys[:stop], rows[:stop]
        hover = ["; ".join(f"{k}: {v}" for k, v in row.items()) for row in rows]
        return _chart(
            key,
            title,
            xs,
            ys,
            assets,
            xlabel=xlabel,
            ylabel=ylabel,
            line=line,
            percent_x=percent_x,
            details=hover,
            csv_name=filename,
        )

    edit_plots = panel(
        "edits-distribution",
        "Edits per fragment",
        editing,
        "edits_per_fragment_csv",
        "edits",
        "fragments",
        xlabel="Edits (last histogram bin may be overflow)",
        ylabel="Fragments",
    )
    edit_plots += panel(
        "rate-distribution",
        "Per-fragment edit rate",
        rates,
        "histogram_csv",
        "bin_start",
        "fragments",
        xlabel="Edit rate (linear ≤1%, logarithmic above)",
        ylabel="Fragments",
        percent_x=True,
    )
    context_items = sorted(metrics.get("context", {}).items())
    context_plot = ""
    if plot:
        context_plot = _chart(
            "context-bias",
            "Trinucleotide context bias",
            list(range(len(context_items))),
            [v.get("edit_fraction") for _, v in context_items],
            assets,
            xlabel="C-centred context",
            ylabel="Edited / opportunities",
            labels=[k for k, _ in context_items],
            details=[
                f"{k}: {v.get('edits')} edits / {v.get('opportunities')} opportunities; fraction {v.get('edit_fraction')}"
                for k, v in context_items
            ],
            csv_name=metrics.get("context_summary", {}).get("csv"),
        )
    if plot and motif_b64 is None and csvs.get(motif.get("pfm_csv")):
        import numpy as np

        from deamtools.qc.qc import _motif_logo_base64

        counts = np.array(
            [[int(float(row[b])) for b in "ACGT"] for row in csvs[motif["pfm_csv"]]]
        )
        counts[len(counts) // 2] = 0
        motif_b64 = _motif_logo_base64(
            counts, motif.get("n_events", 0), motif.get("logo_scale", "frequency")
        )
    logo = "<p class='empty'>Sequence logo not rendered.</p>"
    if plot and motif_b64:
        image_url = assets.image(motif_b64, "sequence-motif")
        logo = f"<img alt='Deaminase motif' src='{image_url}'><div class='downloads'><a class='download' href='{image_url}' download='sequence-motif.png'>↓ PNG</a>{assets.download(motif.get('pfm_csv'))}</div>"
    tss_plot = panel(
        "tss-profile",
        "TSS enrichment profile",
        tss,
        "profile_csv",
        "position",
        "normalized",
        xlabel="Distance from TSS (bp)",
        ylabel="Flank-normalized insertions",
        line=True,
    )
    length_plot = ""
    if layout == "paired-end" and "fragment_length" in metrics:
        length_plot = panel(
            "fragment-length",
            "Fragment length",
            metrics["fragment_length"],
            "histogram_csv",
            "length",
            "pairs",
            xlabel="Fragment length (bp)",
            ylabel="Pairs",
            line=True,
        )

    status_class = {"PASS": "good", "FAIL": "fail"}[assessment["overall"]]
    check_rows = "".join(
        f"<li>{'✓' if passed else '✕'} {_escape(name)}</li>"
        for name, passed in assessment["checks"]
    )
    warnings = "".join(
        f"<li>{_escape(warning)}</li>" for warning in assessment["warnings"]
    )
    overview = (
        f"<section id='overview' class='overall {assessment['overall'].lower()}'><div class='row'><h2>Data checks</h2><span class='badge {status_class}'>{assessment['overall']} {'✓' if assessment['overall']=='PASS' else '✕'}</span></div>"
        f"<div class='row'><strong>{assessment['passed']} / {len(assessment['checks'])} checks passed</strong><span class='badge {'advisory' if assessment['warnings'] else 'observation'}'>Advisories · {len(assessment['warnings'])} notes</span></div><p>{_escape(assessment['summary'])}</p>"
        "<p class='small muted'>Status describes data availability and arithmetic consistency, not assay suitability. No universal editing, duplication or TSS quality thresholds are configured.</p>"
        + _details(
            "Checks and advisory notes", f"<ul>{check_rows}</ul><ul>{warnings}</ul>"
        )
        + "</section>"
    )

    def card(
        label,
        value,
        definition,
        *,
        percent=False,
        compact=False,
        status="Assay-dependent",
        tone="observation",
        suffix="",
    ):
        if not _finite(value):
            status, tone = "N/A", ""
        precise = f"Exact value: {value}. {definition}"
        shown = _number(value, percent=percent, compact=compact) + (
            suffix if _finite(value) else ""
        )
        return (
            f"<div class='card'><div class='card-label'>{label}{_info(precise)}</div>"
            f"<div class='card-value' title='{_escape(precise)}'>{shown}</div>"
            f"<span class='badge {tone}'>{'✓ ' if tone=='good' else '⚠ ' if tone=='advisory' else ''}{status}</span></div>"
        )

    cards = card(
        "Passing reads",
        reads.get("passing"),
        "Sampled records passing flag and MAPQ filters.",
        compact=True,
        status="Available" if reads.get("passing", 0) else "None passing",
        tone="observation" if reads.get("passing", 0) else "fail",
    )
    cards += card(
        "Global edit rate",
        editing.get("global_edit_rate"),
        "Overall fraction of covered editable C/G positions showing C→T or G→A editing. Interpret relative to matched controls and library conditions.",
        percent=True,
    )
    cards += card(
        "Mean edits / fragment",
        editing.get("mean_edits_per_fragment"),
        "Mean count after overlapping mates are merged.",
    )
    cards += card(
        "Mean edit rate / fragment",
        rates.get("mean"),
        "Unweighted mean of per-fragment edit fractions.",
        percent=True,
    )
    cards += card(
        "Duplicate-flagged records",
        reads.get("duplicate_rate"),
        "Duplicate-flagged / analyzed records. This may reflect upstream deduplication and does not measure original library complexity.",
        percent=True,
        status="BAM flag observation",
    )
    score = tss.get("score")
    cards += card(
        "TSS enrichment",
        score,
        "Peak of the flank-normalized profile; flank background = 1× is a normalization reference, not a protocol-independent quality cutoff.",
        status="flank-normalized peak",
        suffix="×",
    )
    overview += f"<div class='cards'>{cards}</div>"

    stats = "".join(
        f"<div><span class='small muted'>{label}</span><strong>{_number(value, percent=percent)}</strong></div>"
        for label, value, percent in (
            ("Mean edits / fragment", editing.get("mean_edits_per_fragment"), False),
            (
                "Median edits / fragment",
                editing.get("median_edits_per_fragment"),
                False,
            ),
            ("Mean edit rate / fragment", rates.get("mean"), True),
            ("Median edit rate / fragment", rates.get("median"), True),
        )
    )
    zero_rows = csvs.get(editing.get("zero_edit_by_opportunities_csv"), [])
    zero_table = "<div class='table-scroll'><table><thead><tr><th>Editable bases</th><th>Fragments</th><th>Zero edits</th><th>Fraction</th></tr></thead><tbody>"
    for row in zero_rows:
        fraction = (
            float(row["zero_edit_fraction"]) if row.get("zero_edit_fraction") else None
        )
        zero_table += f"<tr><td>{_escape(row['editable_bases'])}</td><td>{_number(int(row['fragments']))}</td><td>{_number(int(row['zero_edit_fragments']))}</td><td>{_number(fraction,percent=True)}</td></tr>"
    zero_table += "</tbody></table></div><p class='small muted'>Fragments with no editable bases are separate from fragments with editable bases but no edits. Neither alone establishes closed chromatin.</p>"
    editing_html = (
        "<section id='editing-qc'><h2>Editing signal</h2><p class='subnav'><a href='#editing-statistics'>Editing statistics</a><a href='#fragment-rate'>Per-fragment rate</a><a href='#sequence-motif'>Sequence motif</a><a href='#trinucleotide-context'>Context bias</a></p><div class='grid'><div class='flow'>"
        f"<div><span>Editable C/G</span><div class='number'>{_number(editing.get('total_opportunities'),compact=True)}</div></div><span class='arrow'>→</span>"
        f"<div><span>Edited</span><div class='number'>{_number(editing.get('total_edits'),compact=True)}</div></div><span class='arrow'>→</span>"
        f"<div><span>Global edit rate</span><div class='number'>{_number(editing.get('global_edit_rate'),percent=True)}</div></div></div><div class='mini-grid'>{stats}</div></div>"
        "<p class='small muted'>C→T and G→A are counted regardless of read orientation. Mates are merged; a shared reference position counts once. Rate median is a histogram estimate.</p>"
        f"<div class='grid' id='editing-distributions'>{edit_plots}</div>"
        + _details(
            "Editing statistics", table("editing"), identity="editing-statistics"
        )
        + _details(
            "Per-fragment edit rate",
            table("edit_rate_per_fragment", rate_keys=("mean", "median")),
            identity="fragment-rate",
        )
        + _details("Fragment-level metrics", table("fragments") + zero_table)
        + "</section>"
    )

    bias_stats = {
        "Top context": bias["top"],
        "Top-context edit rate": _number(bias["rate"], percent=True),
        "Overall edit rate": _number(editing.get("global_edit_rate"), percent=True),
        "Context-pooled edit rate": _number(bias["baseline"], percent=True),
        "Top-context / pooled edit rate": (
            _number(bias["fold"]) + "×" if bias["fold"] is not None else "N/A"
        ),
    }
    context_rows = "".join(
        f"<tr><td>{_escape(key)}</td><td>{_number(val.get('edit_fraction'),percent=True)}</td><td>{_number(val.get('edits'))}</td><td>{_number(val.get('opportunities'))}</td></tr>"
        for key, val in context_items
    )
    context_table = (
        "<div class='table-scroll'><table><thead><tr><th>Context</th><th>Edit rate</th><th>Edits</th><th>Opportunities</th></tr></thead><tbody>"
        + context_rows
        + "</tbody></table></div>"
    )
    bias_html = (
        "<section id='enzyme-bias'><h2>Enzyme sequence bias</h2>"
        + _table(
            bias_stats,
            {
                "Top-context / pooled edit rate": "Top-context fraction divided by the opportunity-weighted pooled edit fraction across valid trinucleotide contexts. Descriptive fold difference; not a statistical significance test."
            },
        )
        + "<p class='small muted'>Descriptive fold difference; not a statistical significance test.</p><details><summary>Show detailed motif and context rates</summary><div class='bias-grid'>"
        f"<div id='sequence-motif'><h3>Deaminase sequence motif</h3>{logo}<p class='small muted'>{_number(motif.get('n_events'),compact=True)} contributing events. "
        + (
            "Frequency scale: position 0 is the target cytosine; "
            if motif.get("logo_scale", "frequency") == "frequency"
            else "Bits scale: target centre omitted; "
        )
        + "G→A windows are reverse-complemented.</p></div>"
        f"<div id='trinucleotide-context'>{context_plot}</div></div>"
        "<p class='method'>Top-context fold enrichment uses the opportunity-weighted context-pooled edit rate, so numerator and denominator use the same valid trinucleotide population. It is descriptive, not a statistical significance test. A context with few opportunities can have a high rate. Consider matched controls and bias correction for footprinting; no entropy or KL score is assigned.</p>"
        + _details("Trinucleotide context bias: counts and rates", context_table)
        + _details(
            "Sequence motif metadata",
            table(
                "motif",
                {
                    key: value
                    for key, value in motif.items()
                    if key in {"window", "n_events", "pfm_csv", "logo_scale"}
                },
            ),
        )
        + "</details></section>"
    )

    status = metrics.get("status", {}).get(
        "tss", "not_requested" if not tss else "not_recorded"
    )
    tss_html = (
        "<section id='accessibility-qc'><h2>TSS enrichment</h2>"
        f"<div class='row'><div class='card-value'>{_number(score)}{'×' if _finite(score) else ''}</div><span class='badge {'observation' if _finite(score) else ''}'>{'flank-normalized peak' if _finite(score) else _escape(status.replace('_',' '))}</span></div>"
        "<div class='reference'><strong>Reference: 1× = flank background</strong><br><span class='small'>This is a normalization reference, not an acceptable/ideal cutoff. No assay-validated grading range is configured.</span></div>"
        f"{tss_plot}"
        + _details(
            "TSS methods and measurements",
            f"<p class='method'>Aggregate unshifted read 5′ ends over ±{_plain(tss.get('flank'))} bp around TSS in {_plain(tss.get('bin_size'))}-bp bins; reverse minus-strand profiles and normalize to the two outermost 100-bp flank means (clipped for short windows). The score is the profile maximum. Unlike the ENCODE coverage-based smoothing approach, this implementation counts 1-bp sites directly. Scores depend on annotation, window, bin size and protocol. TSS is not subsampled.</p>"
            + table("tss_enrichment")
            + _table(metrics.get("status", {}).get("tss_annotation", {})),
        )
        + "</section>"
    )

    total, passing = reads.get("total"), reads.get("passing")
    retention = (
        passing / total if _finite(total) and total > 0 and _finite(passing) else None
    )
    filter_rows = "".join(
        f"<tr><td>{label}</td><td title='{_escape(reads[key])}'>{_number(reads[key],compact=True)}</td></tr>"
        for key, label in (
            ("qcfail", "QC-fail flagged"),
            ("low_mapq", "Low MAPQ"),
            ("secondary", "Secondary"),
            ("duplicate", "Duplicate-flagged"),
            ("supplementary", "Supplementary"),
            ("unmapped", "Unmapped"),
        )
        if key in reads
    )
    filter_summary = (
        (
            "<div><h3>Filtering reasons</h3><table><thead><tr><th>Category</th><th>Records</th></tr></thead><tbody>"
            + filter_rows
            + "</tbody></table></div>"
        )
        if filter_rows
        else ""
    )
    alignment = (
        "<section id='alignment-qc'><h2>Read retention</h2><div class='grid'><div><div class='flow'>"
        f"<div>Total analyzed reads<div class='number'>{_number(total,compact=True)}</div></div><span class='arrow'>→</span>"
        f"<div>Passing reads<div class='number'>{_number(passing,compact=True)}</div></div><span class='badge'>{_number(retention,percent=True)} retained</span></div>"
        f"<div class='retention' role='img' aria-label='{_number(retention,percent=True)} reads retained'><span style='width:{max(0,min(100,100*retention)) if retention is not None else 0}%'></span></div>"
        f"</div>{filter_summary}</div><p class='small muted'>Filter categories can overlap and therefore do not sum to the total number of excluded reads. The retention denominator is the same analyzed sample as passing reads, not the unsampled file total.</p>"
        + _details(
            "Read statistics and filtering reasons",
            table("reads"),
            identity="read-statistics",
        )
        + _details("Alignment diagnostics", table("diagnostics"))
        + _details("Full-file counts", table("file_reads"))
    )
    if layout == "paired-end" and "fragment_length" in metrics:
        alignment += _details(
            "Fragment length",
            length_plot + table("fragment_length"),
            identity="fragment-length",
        )
    alignment += "</section>"

    paths = {
        key: provenance[key] for key in ("bam", "fasta", "tss") if provenance.get(key)
    }
    meta = {
        "Sample": sample,
        "Library": layout,
        "BAM": Path(paths["bam"]).name if "bam" in paths else "Not recorded",
        "FASTA": Path(paths["fasta"]).name if "fasta" in paths else "Not recorded",
    }
    if "tss" in paths:
        meta["TSS BED"] = Path(paths["tss"]).name
    meta.update(
        {
            "deamtools": provenance.get("deamtools_version", "Not recorded"),
            "Generated": provenance.get("generated_at", "Not recorded"),
        }
    )
    meta_html = (
        "<table class='meta-table'><tbody>"
        + "".join(
            f"<tr><th>{key}</th><td>{_escape(value)}</td></tr>"
            for key, value in meta.items()
        )
        + "</tbody></table>"
    )
    if layout == "single-end":
        meta_html += "<p class='small muted'>Fragment length: N/A (single-end library). Not applicable to single-end reads.</p>"
    run_html = (
        "<section id='run-information'><h2>Run information</h2>"
        + _details(
            "Parameters and definitions",
            _table(provenance.get("parameters", {}), _PARAMETER_DOCS)
            + _table(
                {
                    "Schema version": metrics.get("schema_version", "legacy"),
                    "Library type": layout,
                    "DeamTools version": provenance.get("deamtools_version"),
                    "Generated": provenance.get("generated_at"),
                    "Editing definition": (
                        "Merged C→T/G→A at covered C/G, independent of read strand"
                        if metrics.get("schema_version") == "2.0"
                        else "Legacy context-restricted global editing counts"
                    ),
                    "TSS window": tss.get("flank"),
                    "TSS bin size": tss.get("bin_size"),
                }
            )
            + table("sampling"),
        )
        + _details(
            "Input filenames (paths redacted)" if redact_paths else "Show full paths",
            _table(paths)
            + (
                "<p class='small muted'>Shareable report: directory components have been removed from provenance and recorded command paths.</p>"
                if redact_paths
                else "<p class='small muted'>Collapsed paths remain stored in this HTML. Use --redact-paths for a shareable report.</p>"
            ),
            identity="full-paths",
        )
    )
    if provenance.get("command"):
        run_html += _details(
            "Command (recorded argv)",
            f"<pre>{_escape(provenance['command'])}</pre><p class='small muted'>Shell-quoted argument vector; original shell whitespace and quoting are not recoverable.</p>",
        )
    run_html += (
        _details(
            "Companion data downloads",
            "<div class='downloads'>"
            + "".join(
                f"<a class='download' href='{href}' download='{_escape(name)}'>{_escape(name)}</a>"
                for name, href in assets.files.items()
            )
            + "</div>",
        )
        + "</section>"
    )
    recommendations = (
        "<section id='recommendations'><h2>QC recommendations</h2><ul class='recommendations'>"
        + "".join(f"<li>{_escape(item)}</li>" for item in assessment["recommendations"])
        + "</ul></section>"
    )
    navigation = [
        ("overview", "Overview"),
        ("editing-qc", "Editing QC"),
        ("editing-distributions", "Distributions"),
        ("enzyme-bias", "Enzyme sequence bias"),
        ("accessibility-qc", "Accessibility QC"),
        ("alignment-qc", "Alignment QC"),
        ("run-information", "Run information"),
        ("recommendations", "Recommendations"),
    ]
    document = (
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{_escape(sample)} · DeamTools QC Report</title><style>{CSS}</style></head><body>"
        "<nav aria-label='Report navigation'><div class='nav-inner'>"
        + "".join(f"<a href='#{key}'>{label}</a>" for key, label in navigation)
        + "</div></nav><main>"
        f"<header><div class='eyebrow'>Chromatin accessibility · Quality control</div><h1>DeamTools QC Report</h1><p class='sample'>{_escape(sample)}</p>"
        + meta_html
        + "</header>"
        + overview
        + editing_html
        + bias_html
        + tss_html
        + alignment
        + run_html
        + recommendations
        + "<noscript><p class='no-js-note'>All measurements and plots are available offline. Enable JavaScript for chart zoom, SVG export and automatic expansion of navigation targets.</p></noscript>"
        + f"<footer class='footer'>Generated by deamtools qc · Data availability is not an assay-quality certification.</footer></main><script>{JS}</script></body></html>"
    )
    target = (
        destination / "report.html"
        if destination
        else Path(html_path) if html_path else source / f"{out_name}.html"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(document, encoding="utf-8")
    return str(target)
