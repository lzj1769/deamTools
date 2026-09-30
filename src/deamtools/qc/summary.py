"""Compare QC reports without pooling distinct libraries or imposing cutoffs."""

from __future__ import annotations

import csv
import html
import json
import math
from pathlib import Path


def run_qc_summary(reports: list[str], out_dir: str, out_name: str) -> list[dict]:
    """Write a CSV and standalone HTML table from versioned QC JSON reports.

    Parameters
    ----------
    reports : list[str]
        Input QC JSON paths. Older reports are retained and labelled legacy.
    out_dir, out_name : str
        Output directory and basename.

    Returns
    -------
    list[dict]
        One comparison row per input report; missing metrics remain null.
    """
    if not reports:
        raise ValueError("At least one QC report is required")
    rows = []
    for filename in reports:
        with open(filename) as stream:
            metrics = json.load(stream)
        if not isinstance(metrics, dict) or not {"editing", "reads"} <= metrics.keys():
            raise ValueError(f"Not a QC report: {filename}")
        provenance = metrics.get("provenance", {})
        parameters = provenance.get("parameters", {})
        tss = metrics.get("tss_enrichment", {})
        row = {
            "sample": provenance.get("sample", Path(filename).stem),
            "report": str(Path(filename).resolve()),
            "schema_version": metrics.get("schema_version", "legacy"),
            "deamtools_version": provenance.get("deamtools_version"),
            "layout": metrics.get("library_layout"),
            "fasta": provenance.get("fasta"),
            "tss_annotation": provenance.get("tss"),
            "min_mapq": parameters.get("min_mapq"),
            "min_baseq": parameters.get("min_baseq"),
            "sampling_fraction": metrics.get("sampling", {}).get("fraction"),
            "passing_reads": metrics["reads"].get("passing"),
            "duplicate_rate": metrics["reads"].get("duplicate_rate"),
            "global_edit_rate": metrics["editing"].get("global_edit_rate"),
            "zero_edit_fraction": metrics["editing"].get("zero_edit_fraction"),
            "mean_edits_per_fragment": metrics["editing"].get(
                "mean_edits_per_fragment"
            ),
            "tss_score": tss.get("score"),
            "tss_flank": tss.get("flank"),
            "tss_status": metrics.get("status", {}).get("tss", "unknown"),
            "non_edit_mismatch_rate": metrics.get("diagnostics", {}).get(
                "non_edit_mismatch_rate"
            ),
        }
        for context, values in sorted(metrics.get("context", {}).items()):
            row[f"context_{context}"] = values.get("edit_fraction")
        row = {
            key: (
                None if isinstance(value, float) and not math.isfinite(value) else value
            )
            for key, value in row.items()
        }
        rows.append(row)
    columns = list(rows[0])
    columns += sorted({key for row in rows for key in row} - set(columns))
    dest = Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / f"{out_name}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    def display(value):
        if value is None or isinstance(value, float) and value != value:
            return "n/a"
        return html.escape(str(value))

    table = (
        "<tr>" + "".join(f"<th>{html.escape(key)}</th>" for key in columns) + "</tr>"
    )
    for row in rows:
        table += (
            "<tr>"
            + "".join(f"<td>{display(row.get(key))}</td>" for key in columns)
            + "</tr>"
        )
    (dest / f"{out_name}.html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>QC comparison</title>"
        "<style>body{font-family:sans-serif;margin:2em}table{border-collapse:collapse}"
        "td,th{border:1px solid #ddd;padding:.5em}th{background:#eee}</style></head>"
        "<body><h1>QC sample comparison</h1><p>Compare matching schema versions, "
        "protocols, references, annotations and filters. Legacy context-restricted "
        "editing rates differ from schema 2.0. Counts can represent different "
        "sampling fractions. No universal pass/fail thresholds are applied.</p>"
        f"<table>{table}</table></body></html>"
    )
    return rows
