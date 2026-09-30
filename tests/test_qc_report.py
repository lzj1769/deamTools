"""Dashboard behavior: interpretation, offline assets, escaping and navigation."""

import json
from html.parser import HTMLParser
from pathlib import Path

import pytest

from deamtools.qc.report import _bias_summary, _number, assess_qc, write_report


class Document(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


@pytest.fixture
def metrics():
    return {
        "schema_version": "2.0",
        "library_layout": "single-end",
        "provenance": {
            "sample": "sample",
            "bam": "/private/input/sample.bam",
            "fasta": "/private/reference/genome.fa",
            "tss": "/private/annotation/tss.bed",
            "generated_at": "2026-09-29T12:00:00+08:00",
            "deamtools_version": "0.1.3",
            "parameters": {"min_mapq": 20, "min_baseq": 20, "threads": 12},
        },
        "reads": {
            "total": 3220000,
            "passing": 2260000,
            "duplicate_rate": 0.03,
            "duplicate": 96600,
            "secondary": 12,
            "qcfail": 34,
            "low_mapq": 56,
        },
        "editing": {
            "total_edits": 21930000,
            "total_opportunities": 191000000,
            "global_edit_rate": 0.1148,
            "mean_edits_per_fragment": 9.703,
            "median_edits_per_fragment": 8,
            "edits_per_fragment_csv": "sample.edits.csv",
        },
        "edit_rate_per_fragment": {
            "mean": 0.12,
            "median": 0.10,
            "histogram_csv": "sample.rate.csv",
        },
        "motif": {"window": 11, "n_events": 42},
        "fragments": {"total": 2260000},
        "context": {
            "ACA": {"edits": 1, "opportunities": 100, "edit_fraction": 0.01},
            "TCA": {"edits": 4, "opportunities": 20, "edit_fraction": 0.2},
        },
        "tss_enrichment": {
            "score": 8.74,
            "n_tss": 40,
            "flank": 2000,
            "bin_size": 10,
            "background": 2,
            "profile_csv": "sample.tss.csv",
        },
        "status": {
            "tss": "ok",
            "tss_annotation": {"input_records": 40, "used_records": 40},
        },
        "sampling": {"subsampled": False, "fraction": 1},
    }


def write_csvs(tmp_path):
    (tmp_path / "sample.edits.csv").write_text(
        "edits,fragments,fraction,is_overflow\n0,10,.5,False\n1,10,.5,False\n"
    )
    (tmp_path / "sample.rate.csv").write_text(
        "bin_start,bin_end,fragments,fraction\n0,.01,4,.4\n.01,.02,6,.6\n"
    )
    (tmp_path / "sample.tss.csv").write_text(
        "position,insertions,mean_insertions_per_tss,normalized\n-5,10,.25,1\n5,20,.5,2\n"
    )


def test_percent_and_human_readable_format():
    assert _number(0.1148, percent=True) == "11.48%"
    assert _number(2260000, compact=True) == "2.26 M"
    assert _number(None) == "N/A"
    assert _number(float("nan")) == "N/A"


def test_descriptive_assessment_does_not_invent_assay_cutoffs(metrics):
    result = assess_qc(metrics)
    assert result["overall"] == "PASS"
    assert result["passed"] == 4
    assert "11.48%" not in result["summary"]
    assert "successfully quantified" in result["summary"]
    assert not any(
        word in result["summary"] for word in ("ideal", "strong", "moderate")
    )
    assert _bias_summary(metrics)["fold"] == pytest.approx(0.2 / (5 / 120))
    metrics["editing"]["global_edit_rate"] = 0.000001
    assert assess_qc(metrics)["passed"] == 4  # not a biological grading threshold
    metrics["reads"]["passing"] = 0
    assert assess_qc(metrics)["overall"] == "FAIL"


def test_tss_missing_background_and_unknown_annotation(metrics):
    metrics["tss_enrichment"]["score"] = None
    metrics["status"]["tss"] = "zero_background"
    metrics["status"]["tss_annotation"]["unknown_contig"] = 2
    result = assess_qc(metrics)
    assert any("zero_background" in item for item in result["warnings"])
    assert any("2 TSS" in item for item in result["warnings"])
    metrics.pop("tss_enrichment")
    metrics["status"] = {"tss": "not_requested"}
    assert not any(
        "TSS enrichment is unavailable" in item
        for item in assess_qc(metrics)["warnings"]
    )


def test_standalone_report_navigation_and_metadata(tmp_path, metrics):
    write_csvs(tmp_path)
    path = write_report(metrics, str(tmp_path), "sample")
    text = Path(path).read_text()
    document = Document(text)
    header = text.split("<header>")[1].split("</header>")[0]
    assert "sample.bam" in header and "/private/" not in header
    assert "11.48%" in text and "2.26 M" in text
    assert "Show full paths" in text and "/private/input/sample.bam" in text
    assert "Fragment length: N/A (single-end library)" in text
    assert "id='fragment-length'" not in text
    assert "other_filtered" not in text
    assert "<th>Info</th>" in text and "<th>Meaning</th>" not in text
    assert "data:text/csv;base64," in text
    ids = {attrs["id"] for _, attrs in document.tags if "id" in attrs}
    for tag, attrs in document.tags:
        if tag == "a" and attrs.get("href", "").startswith("#"):
            assert attrs["href"][1:] in ids
        if tag in ("script", "img", "link"):
            assert not attrs.get("src", attrs.get("href", "")).startswith(
                ("http:", "https:")
            )
    assert {"editing-qc", "enzyme-bias", "alignment-qc", "run-information"} <= ids
    assert "data-chart='edits-distribution'" in text
    assert "data-chart='rate-distribution'" in text
    assert "data-chart='tss-profile'" in text


def test_directory_report_has_portable_assets_and_csvs(tmp_path, metrics):
    write_csvs(tmp_path)
    dest = tmp_path / "report"
    # A minimal payload is enough to test externalization; logo rendering is tested separately.
    path = write_report(
        metrics, str(tmp_path), "sample", motif_b64="aGVsbG8=", report_dir=str(dest)
    )
    assert path == str(dest / "report.html")
    text = Path(path).read_text()
    assert "data:image/png;base64," not in text
    assert "data:text/csv;base64," not in text
    assert (dest / "assets/sequence-motif.png").read_bytes() == b"hello"
    assert (dest / "assets/tss-profile.svg").is_file()
    assert (dest / "sample.tss.csv").read_bytes() == (
        tmp_path / "sample.tss.csv"
    ).read_bytes()
    assert json.loads((dest / "metrics.json").read_text()) == metrics
    for tag, attrs in Document(text).tags:
        if tag == "img":
            assert (dest / attrs["src"]).is_file()


def test_paths_commands_and_sample_names_are_escaped(tmp_path, metrics):
    attack = "<script>alert('x')</script>"
    metrics["provenance"]["sample"] = attack
    metrics["provenance"]["command"] = f"deamtools qc --out_name {attack}"
    metrics["provenance"]["bam"] = "/data/<unsafe>.bam"
    text = Path(write_report(metrics, str(tmp_path), "sample", plot=False)).read_text()
    assert attack not in text
    assert "&lt;script&gt;" in text and "&lt;unsafe&gt;.bam" in text
    assert "Command (recorded argv)" in text
    assert "<svg" not in text and "data:image/" not in text
    assert len([tag for tag, _ in Document(text).tags if tag == "script"]) == 1


def test_paired_length_and_missing_legacy_provenance(tmp_path, metrics):
    metrics["library_layout"] = "paired-end"
    metrics["fragment_length"] = {"mean": 210, "median": 180, "n_pairs": 32}
    metrics.pop("provenance")
    text = Path(write_report(metrics, str(tmp_path), "sample", plot=False)).read_text()
    assert "id='fragment-length'" in text
    assert "Command (recorded argv)" not in text
    assert "Not recorded" in text
    assert "Fragment length: N/A" not in text


def test_report_does_not_read_companion_paths_outside_source(tmp_path, metrics):
    metrics["editing"]["edits_per_fragment_csv"] = "../private.csv"
    text = Path(write_report(metrics, str(tmp_path), "sample", plot=False)).read_text()
    assert "download='../private.csv'" not in text


def test_integrity_pass_separate_from_advisories_and_observations(tmp_path, metrics):
    metrics["status"]["tss_annotation"]["unknown_contig"] = 2
    metrics["reads"]["duplicate_rate"] = 0
    metrics["reads"]["duplicate"] = 0
    text = Path(write_report(metrics, str(tmp_path), "sample", plot=False)).read_text()
    assert "PASS ✓" in text and "4 / 4 checks passed" in text
    assert "Advisories · 2 notes" in text
    assert "WARNING · data checks" not in text
    assert text.count("class='badge good'") == 1  # validation only
    assert "No duplicates observed" not in text
    assert "No duplicate-flagged records were present" in text
    assert "This may reflect upstream deduplication." in text
    assert "Duplicate-flagged records" in text
    assert "8.74×" in text and "flank-normalized peak" in text
    assert "Above flank background" not in text
    summary = assess_qc(metrics)["summary"]
    assert "11.48" not in summary and "8.74" not in summary


def test_bias_summary_and_filter_reasons_are_visible_by_default(tmp_path, metrics):
    text = Path(write_report(metrics, str(tmp_path), "sample", plot=False)).read_text()
    bias = text.split("id='enzyme-bias'", 1)[1].split("</section>", 1)[0]
    visible = bias.split("<details>", 1)[0]
    assert "Top context" in visible and "Top-context / pooled edit rate" in visible
    assert (
        "Descriptive fold difference; not a statistical significance test." in visible
    )
    assert "Show detailed motif and context rates" in bias
    assert "Top / context-pooled rate" not in text
    alignment = text.split("id='alignment-qc'", 1)[1].split("<details", 1)[0]
    for label in ("QC-fail flagged", "Low MAPQ", "Secondary", "Duplicate-flagged"):
        assert label in alignment
    assert "Filter categories can overlap and therefore do not sum" in alignment
    assert "Other filtered" not in alignment


def test_tooltips_add_information_or_are_absent(tmp_path, metrics):
    from deamtools.qc.qc import _METRIC_DOCS
    from deamtools.qc.report import _table

    text = Path(
        write_report(
            metrics, str(tmp_path), "sample", plot=False, descriptions=_METRIC_DOCS
        )
    ).read_text()
    assert "Minimum MAPQ required for a mapped read" in text
    assert "Minimum base quality required when counting editable C/G" in text
    assert "single most important" not in text
    assert "Interpret relative to matched controls and library conditions" in text
    assert "class='info'" not in _table({"bam": "a.bam", "unexplained": 3})


@pytest.mark.parametrize("directory_mode", [False, True])
def test_shareable_report_removes_paths_without_mutating_original(
    tmp_path, metrics, directory_mode
):
    import copy
    import shlex

    secret = "/private/research user/project"
    metrics["provenance"].update(
        {
            "bam": secret + "/sample.bam",
            "fasta": secret + "/genome.fa",
            "tss": secret + "/tss.bed",
        }
    )
    metrics["provenance"]["parameters"]["report_dir"] = secret + "/dashboard"
    metrics["provenance"]["command"] = shlex.join(
        [
            secret + "/.venv/bin/deamtools",
            "qc",
            "--bam",
            secret + "/sample.bam",
            "--fasta=" + secret + "/genome.fa",
            "--tss",
            secret + "/tss.bed",
            "--out_dir",
            secret + "/output",
            "--report-dir",
            secret + "/dashboard",
        ]
    )
    original = copy.deepcopy(metrics)
    destination = tmp_path / "shareable" if directory_mode else None
    path = write_report(
        metrics,
        str(tmp_path),
        "sample",
        plot=False,
        report_dir=str(destination) if destination else None,
        redact_paths=True,
    )
    text = Path(path).read_text()
    assert secret not in text
    assert "/private/" not in text
    assert "deamtools qc --bam sample.bam" in text
    assert "--fasta=genome.fa" in text
    assert "Input filenames (paths redacted)" in text
    assert "sample.bam" in text
    assert metrics == original
    if destination:
        copy_text = (destination / "metrics.json").read_text()
        assert secret not in copy_text
        assert json.loads(copy_text)["editing"] == original["editing"]
    full_path = write_report(metrics, str(tmp_path), "full", plot=False)
    assert secret in Path(full_path).read_text()


def test_redaction_handles_unparseable_and_windows_commands():
    from deamtools.qc.report import _redacted_command

    assert "/secret" not in _redacted_command("deamtools --bam '/secret/unclosed")
    assert (
        _redacted_command(r"deamtools --bam 'C:\private\sample.bam'")
        == "deamtools --bam sample.bam"
    )
