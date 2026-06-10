"""Tests for deamtools.stat.qc."""

import json
import os

import pysam
import pytest

from deamtools.qc import run_qc

# Reference: A C G T C G A T C G   (positions 0..9)
# Forward C positions: 1, 4, 8   Reverse G positions: 2, 5, 9
REF_SEQ = "ACGTCGATCG"

BAM_HEADER = {
    "HD": {"VN": "1.6"},
    "SQ": [{"LN": len(REF_SEQ), "SN": "chr1"}],
}


def _make_read(
    name: str,
    seq: str,
    pos: int,
    *,
    is_reverse: bool = False,
    is_paired: bool = True,
    is_read1: bool = True,
    mate_pos: int | None = None,
    mate_reverse: bool = False,
    proper_pair: bool = True,
    template_length: int = 0,
    mapq: int = 30,
    baseq: int = 40,
    extra_flags: int = 0,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    flag = extra_flags
    if is_reverse:
        flag |= 0x10
    if is_paired:
        flag |= 0x1
        if proper_pair:
            flag |= 0x2
        if mate_reverse:
            flag |= 0x20
        flag |= 0x40 if is_read1 else 0x80
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * len(seq))
    a.template_length = template_length
    if is_paired and mate_pos is not None:
        a.next_reference_id = 0
        a.next_reference_start = mate_pos
    return a


def _write_bam(path: str, reads: list) -> str:
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=BAM_HEADER) as bam:
        for read in reads:
            bam.write(read)
    pysam.sort("-o", path, tmp)
    os.remove(tmp)
    pysam.index(path)
    return path


@pytest.fixture()
def fasta_file(tmp_path):
    path = str(tmp_path / "ref.fa")
    with open(path, "w") as f:
        f.write(f">chr1\n{REF_SEQ}\n")
    pysam.faidx(path)
    return path


class TestQC:
    def test_edit_rate_and_opportunities(self, tmp_path, fasta_file):
        # One forward read covering pos 0-9 with a single C->T edit at pos 4.
        # Edit calling is strand-agnostic (matching bam2bw): every reference C
        # or G at a position with both flanks present (internal pos 1..8) is an
        # opportunity. Ref ACGTCGATCG -> C/G at 1,2,4,5,8 = 5 opportunities.
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False)  # T at pos 4
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)

        assert m["editing"]["total_opportunities"] == 5  # C/G at 1,2,4,5,8
        assert m["editing"]["total_edits"] == 1
        assert m["editing"]["global_edit_rate"] == pytest.approx(1 / 5)

    def test_reverse_read_g_to_a_counted(self, tmp_path, fasta_file):
        # Reverse read covering pos 1-9 with G->A at pos 5.
        # Strand-agnostic: opportunities are reference C or G at internal
        # positions 1..8 -> 1(C),2(G),4(C),5(G),8(C) = 5; one G->A edit at pos 5.
        # ref[1:10] = CGTCGATCG ; change pos5 G->A: CGTCAATCG
        read = _make_read(
            "r1", "CGTCAATCG", 1, is_reverse=True, is_paired=False
        )
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)

        assert m["editing"]["total_opportunities"] == 5  # C/G at 1,2,4,5,8
        assert m["editing"]["total_edits"] == 1

    def test_forward_read_g_to_a_counted(self, tmp_path, fasta_file):
        # Strand-agnostic: a G->A mismatch on a FORWARD read is now an edit
        # (previously, forward reads only counted C->T). ref pos 2 = G -> A.
        # ref ACGTCGATCG -> ACATCGATCG.
        read = _make_read("r1", "ACATCGATCG", 0, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        assert m["editing"]["total_edits"] == 1
        # The G->A edit's context is reverse-complemented to C-centered: ref[1:4]
        # = "CGT" -> revcomp "ACG".
        assert m["context"]["ACG"]["edits"] == 1

    def test_context_is_c_centered_for_both_strands(self, tmp_path, fasta_file):
        # Forward C->T edit at pos 4: ref context = ref[3:6] = "TCG".
        fwd = _make_read("f", "ACGTTGATCG", 0, is_paired=False)  # T at pos 4
        # Reverse G->A edit at pos 5: ref[4:7]="CGA", revcomp="TCG".
        # ref[1:10]=CGTCGATCG, change pos5 G->A -> CGTCAATCG
        rev = _make_read("r", "CGTCAATCG", 1, is_reverse=True, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [fwd, rev])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)

        # Both edits land in the unified C-centered "TCG" context.
        assert "TCG" in m["context"]
        assert m["context"]["TCG"]["edits"] == 2

    def test_edit_rate_per_read(self, tmp_path, fasta_file):
        # Forward read over the whole reference ACGTCGATCG.
        # Editable bases (C or G), strand-agnostic: C@1, G@2, C@4, G@5, C@8, G@9
        #   -> 6 editable bases.
        # Read "ACGTTGATCG" edits C->T at pos 4 only -> 1 edited.
        # Per-read edit rate = 1/6 ~= 0.1667.
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)

        erpr = m["edit_rate_per_read"]
        assert erpr["n_reads"] == 1
        assert erpr["mean"] == pytest.approx(1 / 6, abs=1e-6)
        assert sum(erpr["histogram"]) == 1
        assert len(erpr["histogram"]) == len(erpr["bin_edges"]) - 1
        # The single read falls in the bin covering 1/6.
        assert erpr["histogram"][int((1 / 6) * len(erpr["histogram"]))] == 1

    def test_edit_rate_counts_both_strands_as_editable(self, tmp_path, fasta_file):
        # A reverse read's editable bases are also counted as reference C or G.
        # ref[1:10] = CGTCGATCG -> C/G at every position except T@3, A@6, T@7
        #   editable = C@1,G@2,C@4,G@5,C@8,G@9 = 6; one G->A edit at pos 5.
        rev = _make_read("r", "CGTCAATCG", 1, is_reverse=True, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [rev])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        assert m["edit_rate_per_read"]["mean"] == pytest.approx(1 / 6, abs=1e-6)

    def test_min_baseq_excludes_opportunity(self, tmp_path, fasta_file):
        # Low base quality everywhere -> no opportunities counted at all.
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False, baseq=5)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=20, plot=False)
        assert m["editing"]["total_opportunities"] == 0
        assert m["editing"]["total_edits"] == 0

    def test_read_counts_and_duplicate_rate(self, tmp_path, fasta_file):
        good = _make_read("r1", REF_SEQ, 0, is_paired=False)
        dup = _make_read("r2", REF_SEQ, 0, is_paired=False, extra_flags=0x400)
        bam = _write_bam(str(tmp_path / "x.bam"), [good, dup])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        assert m["reads"]["total"] == 2
        assert m["reads"]["duplicate"] == 1
        assert m["reads"]["duplicate_rate"] == pytest.approx(0.5)
        assert m["reads"]["passing"] == 1  # duplicate is filtered out

    def test_min_mapq_filters_read(self, tmp_path, fasta_file):
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False, mapq=10)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=20, min_baseq=0, plot=False)
        assert m["reads"]["passing"] == 0
        assert m["editing"]["total_opportunities"] == 0

    def test_fragment_length_from_proper_pair_read1(self, tmp_path, fasta_file):
        r1 = _make_read(
            "p1", REF_SEQ[0:6], 0, is_read1=True, mate_pos=4,
            mate_reverse=True, template_length=10,
        )
        r2 = _make_read(
            "p1", REF_SEQ[4:10], 4, is_read1=False, is_reverse=True,
            mate_pos=0, template_length=-10,
        )
        bam = _write_bam(str(tmp_path / "x.bam"), [r1, r2])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        # Only read1 contributes one fragment of length 10.
        assert m["fragment_length"]["n_pairs"] == 1
        assert m["fragment_length"]["median"] == pytest.approx(10)

    def test_json_and_html_written(self, tmp_path, fasta_file):
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path / "sub")  # nested dir is created
        run_qc(bam, fasta_file, out_dir, "sample", min_mapq=0, min_baseq=0, plot=True)
        json_path = os.path.join(out_dir, "sample.json")
        html_path = os.path.join(out_dir, "sample.html")
        assert os.path.exists(json_path)
        assert os.path.exists(html_path)
        with open(json_path) as f:
            data = json.load(f)
        assert "editing" in data and "context" in data
        html = open(html_path).read()
        assert "DeamTools QC Report" in html
        # Plot embedded and metric descriptions present.
        assert "data:image/png;base64," in html
        assert "global_edit_rate" in html and "Meaning" in html

    def test_html_omits_image_when_no_plot(self, tmp_path, fasta_file):
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        html = open(os.path.join(out_dir, "qc.html")).read()
        assert "data:image/png;base64," not in html
        # Tables and descriptions are still present without the figure.
        assert "Trinucleotide context bias" in html

    def test_tss_enrichment_computed(self, tmp_path, fasta_file):
        # Pile reads so insertion 5' ends concentrate at the TSS center (pos 5).
        reads = [
            _make_read(f"r{i}", REF_SEQ[5:10], 5, is_paired=False)
            for i in range(20)
        ]
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        tss = str(tmp_path / "tss.bed")
        with open(tss, "w") as f:
            f.write("chr1\t4\t6\n")  # midpoint = 5
        out_dir = str(tmp_path)
        m = run_qc(
            bam, fasta_file, out_dir, "qc", tss_path=tss, min_mapq=0,
            min_baseq=0, tss_flank=4, plot=False,
        )
        assert "tss_enrichment" in m
        assert len(m["tss_enrichment"]["profile"]) == 2 * 4 + 1
