"""Tests for deamtools.qc."""

import csv
import json
import os

import numpy as np
import pysam
import pytest

from deamtools.qc import qc, run_qc

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
        read = _make_read("r1", "CGTCAATCG", 1, is_reverse=True, is_paired=False)
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

    def test_edit_rate_per_fragment(self, tmp_path, fasta_file):
        # Forward read over the whole reference ACGTCGATCG.
        # Editable bases (C or G), strand-agnostic: C@1, G@2, C@4, G@5, C@8, G@9
        #   -> 6 editable bases.
        # Read "ACGTTGATCG" edits C->T at pos 4 only -> 1 edited.
        # Per-read edit rate = 1/6 ~= 0.1667.
        read = _make_read("r1", "ACGTTGATCG", 0, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)

        erpr = m["edit_rate_per_fragment"]
        assert erpr["n_fragments_with_editable_bases"] == 1
        assert erpr["mean"] == pytest.approx(1 / 6, abs=1e-6)
        # The histogram lives in its CSV, like the TSS profile.
        assert "histogram" not in erpr
        with open(os.path.join(out_dir, erpr["histogram_csv"])) as f:
            rows = list(csv.DictReader(f))
        assert sum(int(r["fragments"]) for r in rows) == 1
        # The single fragment falls in the bin covering 1/6.
        hit = [r for r in rows if int(r["fragments"])]
        assert float(hit[0]["bin_start"]) <= 1 / 6 < float(hit[0]["bin_end"])

    def test_edit_rate_counts_both_strands_as_editable(self, tmp_path, fasta_file):
        # A reverse read's editable bases are also counted as reference C or G.
        # ref[1:10] = CGTCGATCG -> C/G at every position except T@3, A@6, T@7
        #   editable = C@1,G@2,C@4,G@5,C@8,G@9 = 6; one G->A edit at pos 5.
        rev = _make_read("r", "CGTCAATCG", 1, is_reverse=True, is_paired=False)
        bam = _write_bam(str(tmp_path / "x.bam"), [rev])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        assert m["edit_rate_per_fragment"]["mean"] == pytest.approx(1 / 6, abs=1e-6)

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
            "p1",
            REF_SEQ[0:6],
            0,
            is_read1=True,
            mate_pos=4,
            mate_reverse=True,
            template_length=10,
        )
        r2 = _make_read(
            "p1",
            REF_SEQ[4:10],
            4,
            is_read1=False,
            is_reverse=True,
            mate_pos=0,
            template_length=-10,
        )
        bam = _write_bam(str(tmp_path / "x.bam"), [r1, r2])
        out_dir = str(tmp_path)
        m = run_qc(bam, fasta_file, out_dir, "qc", min_mapq=0, min_baseq=0, plot=False)
        # Only read1 contributes one fragment of length 10.
        assert m["fragment_length"]["n_pairs"] == 1
        assert m["fragment_length"]["median"] == pytest.approx(10)

    def test_motif_logo_built_from_bam(self, tmp_path):
        # Reference long enough for the 11-bp motif window around a central C.
        ref = "A" * 20 + "C" + "A" * 20  # C at index 20
        refpath = str(tmp_path / "mref.fa")
        with open(refpath, "w") as f:
            f.write(f">chr1\n{ref}\n")
        pysam.faidx(refpath)

        read_seq = "A" * 20 + "T" + "A" * 20  # C->T edit at index 20
        hdr = {"HD": {"VN": "1.6"}, "SQ": [{"LN": len(ref), "SN": "chr1"}]}
        a = pysam.AlignedSegment()
        a.query_name = "r1"
        a.query_sequence = read_seq
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 30
        a.cigar = [(0, len(read_seq))]
        a.query_qualities = pysam.qualitystring_to_array("I" * len(read_seq))
        unsorted = str(tmp_path / "m.unsorted.bam")
        bam = str(tmp_path / "m.bam")
        with pysam.AlignmentFile(unsorted, "wb", header=hdr) as out:
            out.write(a)
        pysam.sort("-o", bam, unsorted)
        pysam.index(bam)

        out_dir = str(tmp_path)
        m = run_qc(bam, refpath, out_dir, "mqc", min_mapq=0, min_baseq=0, plot=True)
        assert m["motif"]["window"] == 11
        assert m["motif"]["n_events"] == 1  # one C->T edit with a full window
        html = open(os.path.join(out_dir, "mqc.html")).read()
        assert "Deaminase sequence motif" in html
        # Two embedded PNGs: the summary panel figure and the motif logo.
        assert html.count("data:image/png;base64,") >= 2

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


class TestDistributionCsvs:
    """The numbers behind the summary figure's histograms, one CSV each."""

    def _run(self, tmp_path, fasta_file, reads, **kw):
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "o")
        m = run_qc(bam, fasta_file, out, "s", min_mapq=0, min_baseq=0, **kw)
        return m, out

    def _rows(self, out, name):
        with open(os.path.join(out, name)) as f:
            return list(csv.DictReader(f))

    def test_edits_per_fragment_csv(self, tmp_path, fasta_file):
        # Two fragments with one C->T edit, one with none.
        reads = [
            _make_read("a", "ACGTTGATCG", 0, is_paired=False),
            _make_read("b", "ACGTTGATCG", 0, is_paired=False),
            _make_read("c", REF_SEQ, 0, is_paired=False),
        ]
        m, out = self._run(tmp_path, fasta_file, reads, plot=False)
        name = m["editing"]["edits_per_fragment_csv"]
        assert name == "s.edits_per_fragment.csv"
        rows = self._rows(out, name)
        assert list(rows[0]) == ["edits", "fragments", "fraction", "is_overflow"]
        assert len(rows) == qc._MAX_EDITS + 1
        assert int(rows[0]["fragments"]) == 1 and int(rows[1]["fragments"]) == 2
        assert float(rows[1]["fraction"]) == pytest.approx(2 / 3, rel=1e-5)
        # Only the last row is the overflow bin.
        assert [r["is_overflow"] for r in rows].count("True") == 1
        assert rows[-1]["is_overflow"] == "True"

    def test_edit_rate_csv_bins_tile_zero_to_one(self, tmp_path, fasta_file):
        reads = [_make_read("a", "ACGTTGATCG", 0, is_paired=False)]
        m, out = self._run(tmp_path, fasta_file, reads, plot=False)
        name = m["edit_rate_per_fragment"]["histogram_csv"]
        assert name == "s.edit_rate_per_fragment.csv"
        rows = self._rows(out, name)
        assert list(rows[0]) == ["bin_start", "bin_end", "fragments", "fraction"]
        assert len(rows) == qc._RATE_BINS
        assert float(rows[0]["bin_start"]) == 0.0
        assert float(rows[-1]["bin_end"]) == 1.0
        for prev, cur in zip(rows, rows[1:], strict=False):
            assert float(prev["bin_end"]) == pytest.approx(float(cur["bin_start"]))

    def test_csvs_are_written_even_without_plots(self, tmp_path, fasta_file):
        reads = [_make_read("a", REF_SEQ, 0, is_paired=False)]
        _, out = self._run(tmp_path, fasta_file, reads, plot=False)
        assert os.path.exists(os.path.join(out, "s.edits_per_fragment.csv"))
        assert os.path.exists(os.path.join(out, "s.edit_rate_per_fragment.csv"))
        # No --tss, so no TSS CSV.
        assert not os.path.exists(os.path.join(out, "s.tss_enrichment.csv"))

    def test_report_header_is_a_table(self, tmp_path, fasta_file):
        reads = [_make_read("a", REF_SEQ, 0, is_paired=False)]
        _, out = self._run(tmp_path, fasta_file, reads, plot=False)
        html = open(os.path.join(out, "s.html")).read()
        assert "<table class='meta-table'>" in html
        for label in ("Sample", "Library", "BAM", "FASTA", "deamtools", "Generated"):
            assert f"<tr><th>{label}</th>" in html
        assert "<th>TSS BED</th>" not in html  # only listed when --tss is given


class TestMotif:
    """The motif's counts go to CSV; the logo's y axis is bits or frequency."""

    # Two C->T edits at pos 4, whose 11-bp window runs past both ends of the
    # 10-bp reference -- so use a longer reference for these.
    REF = "AATTCCGGAATTCCGGAATTCCGGAATT"

    def _setup(self, tmp_path):
        fasta = str(tmp_path / "long.fa")
        with open(fasta, "w") as f:
            f.write(f">chr1\n{self.REF}\n")
        pysam.faidx(fasta)
        header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": len(self.REF), "SN": "chr1"}]}
        edited = self.REF[:12] + "T" + self.REF[13:]  # ref[12] is C -> T
        path = str(tmp_path / "m.bam")
        tmp = path + ".u.bam"
        with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
            for i in range(2):
                bam.write(_make_read(f"r{i}", edited, 0, is_paired=False))
        pysam.sort("-o", path, tmp)
        pysam.index(path)
        return path, fasta

    def test_pfm_csv_holds_counts_with_the_target_c(self, tmp_path):
        bam, fasta = self._setup(tmp_path)
        out = str(tmp_path / "o")
        m = run_qc(bam, fasta, out, "s", min_mapq=0, min_baseq=0, plot=False)
        assert m["motif"]["n_events"] == 2
        assert m["motif"]["pfm_csv"] == "s.motif_pfm.csv"
        with open(os.path.join(out, m["motif"]["pfm_csv"])) as f:
            rows = {int(r["position"]): r for r in csv.DictReader(f)}
        assert sorted(rows) == list(range(-5, 6))
        # Position 0 is the edited base: always C in the unified orientation.
        assert [int(rows[0][b]) for b in "ACGT"] == [0, 2, 0, 0]
        # Each flank row puts both events on the reference base at that offset
        # from the edit (ref[12]); derived from REF rather than hand-counted.
        for d in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5):
            assert int(rows[d][self.REF[12 + d]]) == 2, d
        for pos, r in rows.items():
            assert sum(int(r[b]) for b in "ACGT") == 2, pos

    def test_both_scales_render(self, tmp_path):
        bam, fasta = self._setup(tmp_path)
        for scale in ("bits", "frequency"):
            out = str(tmp_path / scale)
            m = run_qc(
                bam,
                fasta,
                out,
                "s",
                min_mapq=0,
                min_baseq=0,
                plot=True,
                logo_scale=scale,
            )
            assert m["motif"]["logo_scale"] == scale
            html = open(os.path.join(out, "s.html")).read()
            assert "alt='Deaminase motif'" in html
            assert "s.motif_pfm.csv" in html
        freq_html = open(os.path.join(tmp_path / "frequency", "s.html")).read()
        assert "position 0 is the target cytosine" in freq_html

    def test_frequency_matrix_rows_sum_to_one(self):
        pwm = np.array([[3, 1, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1]])
        counts = qc._motif_counts_df(pwm, n_events=4)
        assert list(counts.index) == [-1, 0, 1]
        assert list(counts.loc[0]) == [0, 4, 0, 0]
        freq = counts.div(counts.sum(axis=1), axis=0)
        assert np.allclose(freq.sum(axis=1), 1.0)

    def test_rejects_an_unknown_scale(self, tmp_path):
        bam, fasta = self._setup(tmp_path)
        with pytest.raises(ValueError, match="logo_scale"):
            run_qc(bam, fasta, str(tmp_path / "o"), "s", logo_scale="percent")


class TestLibraryLayout:
    """Pair-dependent metrics are reported only for a paired-end library.

    A single-end BAM has no insert size and no pairs, so a proper-pair rate of 0
    or a fragment length of 0 bp would read as a failure rather than as "does
    not apply". The layout is decided from the head of the file before the main
    pass.
    """

    EDITED = "ACGTTGATCG"  # one C->T at pos 4

    def _single(self, tmp_path, n=3):
        reads = [_make_read(f"r{i}", self.EDITED, 0, is_paired=False) for i in range(n)]
        return _write_bam(str(tmp_path / "se.bam"), reads)

    def _paired(self, tmp_path):
        reads = [
            _make_read(
                "p",
                self.EDITED,
                0,
                is_read1=True,
                mate_pos=0,
                mate_reverse=True,
                template_length=10,
            ),
            _make_read(
                "p",
                self.EDITED,
                0,
                is_read1=False,
                is_reverse=True,
                mate_pos=0,
                template_length=-10,
            ),
        ]
        return _write_bam(str(tmp_path / "pe.bam"), reads)

    def test_detects_single_end(self, tmp_path):
        assert qc._detect_layout(self._single(tmp_path)) == qc.LAYOUT_SINGLE

    def test_detects_paired_end(self, tmp_path):
        assert qc._detect_layout(self._paired(tmp_path)) == qc.LAYOUT_PAIRED

    def test_one_paired_record_makes_the_library_paired(self, tmp_path):
        reads = [_make_read(f"s{i}", REF_SEQ, 0, is_paired=False) for i in range(5)]
        reads.append(_make_read("p", REF_SEQ, 0, mate_pos=0))
        bam = _write_bam(str(tmp_path / "mixed.bam"), reads)
        assert qc._detect_layout(bam) == qc.LAYOUT_PAIRED

    def test_empty_bam_is_single_end(self, tmp_path):
        bam = _write_bam(str(tmp_path / "empty.bam"), [])
        assert qc._detect_layout(bam) == qc.LAYOUT_SINGLE

    def test_single_end_omits_pair_metrics(self, tmp_path, fasta_file):
        m = run_qc(
            self._single(tmp_path),
            fasta_file,
            str(tmp_path / "o"),
            "se",
            min_mapq=0,
            min_baseq=0,
            plot=False,
        )
        assert m["library_layout"] == "single-end"
        assert "proper_pair" not in m["reads"]
        assert "proper_pair_rate" not in m["reads"]
        assert "fragment_length" not in m
        # Everything that does apply is still there.
        assert m["reads"]["passing"] == 3
        assert m["fragments"]["total"] == 3
        assert m["editing"]["total_edits"] == 3

    def test_paired_end_keeps_pair_metrics(self, tmp_path, fasta_file):
        m = run_qc(
            self._paired(tmp_path),
            fasta_file,
            str(tmp_path / "o"),
            "pe",
            min_mapq=0,
            min_baseq=0,
            plot=False,
        )
        assert m["library_layout"] == "paired-end"
        assert m["reads"]["proper_pair"] == 2
        assert m["reads"]["proper_pair_rate"] == pytest.approx(1.0)
        assert m["fragment_length"]["n_pairs"] == 1

    def test_single_end_report_says_not_applicable(self, tmp_path, fasta_file):
        out = str(tmp_path / "o")
        run_qc(
            self._single(tmp_path),
            fasta_file,
            out,
            "se",
            min_mapq=0,
            min_baseq=0,
            plot=True,
        )
        html = open(os.path.join(out, "se.html")).read()
        assert "single-end" in html
        assert "proper_pair_rate" not in html
        assert "Not applicable" in html
        # The JSON on disk matches what run_qc returned.
        on_disk = json.load(open(os.path.join(out, "se.json")))
        assert on_disk["library_layout"] == "single-end"
        assert "fragment_length" not in on_disk


class TestFragmentMerging:
    """Editing is counted per fragment: overlapping mates are one observation.

    REF_SEQ is ACGTCGATCG. Editable C/G sit at 1, 2, 4, 5, 8, 9; the five with
    both flanking bases present (1, 2, 4, 5, 8) are context opportunities.
    """

    # Both mates of a real pair carry the same edits: library prep turns the
    # deaminated C into a T:A pair that both strands then report.
    EDITED = "ATGTCAATCG"  # C->T at pos 1, G->A at pos 5

    def _pair(self, seq1, seq2, **kw):
        return [
            _make_read(
                "frag", seq1, 0, is_read1=True, mate_pos=0, mate_reverse=True, **kw
            ),
            _make_read(
                "frag", seq2, 0, is_read1=False, is_reverse=True, mate_pos=0, **kw
            ),
        ]

    def _run(self, tmp_path, fasta_file, reads, name="f", **kw):
        bam = _write_bam(str(tmp_path / f"{name}.bam"), reads)
        return run_qc(
            bam,
            fasta_file,
            str(tmp_path / f"o_{name}"),
            name,
            min_baseq=0,
            plot=False,
            **kw,
        )

    def test_overlapping_mates_count_a_position_once(self, tmp_path, fasta_file):
        m = self._run(
            tmp_path, fasta_file, self._pair(self.EDITED, self.EDITED), min_mapq=0
        )
        assert m["reads"]["total"] == 2
        assert m["fragments"]["total"] == 1
        assert m["fragments"]["from_mate_pairs"] == 1

        # Per record this was 10 opportunities and 4 edits; the mates cover the
        # same 0-9, so the fragment sees each position exactly once.
        assert m["editing"]["total_opportunities"] == 5
        assert m["editing"]["total_edits"] == 2
        assert m["editing"]["mean_edits_per_fragment"] == pytest.approx(2.0)
        # 2 edited of 6 editable C/G, not 4 of 12.
        assert m["edit_rate_per_fragment"]["mean"] == pytest.approx(2 / 6, abs=1e-6)
        assert m["edit_rate_per_fragment"]["n_fragments_with_editable_bases"] == 1

    def test_non_overlapping_mates_still_add_up(self, tmp_path, fasta_file):
        # R1 covers 0-4, R2 covers 5-9: disjoint, so nothing is deduplicated and
        # the fragment sees the union.
        reads = [
            _make_read(
                "frag", "ATGTC", 0, is_read1=True, mate_pos=5, mate_reverse=True
            ),  # ref[0:5] ACGTC, C->T at 1
            _make_read(
                "frag", "AATCG", 5, is_read1=False, is_reverse=True, mate_pos=0
            ),  # ref[5:10] GATCG, G->A at 5
        ]
        m = self._run(tmp_path, fasta_file, reads, min_mapq=0)
        assert m["fragments"]["total"] == 1
        assert m["editing"]["total_opportunities"] == 5  # 1,2,4 from R1; 5,8 from R2
        assert m["editing"]["total_edits"] == 2
        assert m["edit_rate_per_fragment"]["mean"] == pytest.approx(2 / 6, abs=1e-6)

    def test_higher_base_quality_wins_a_disagreement(self, tmp_path, fasta_file):
        """At an overlap the mates disagreeing means one of them misread."""

        def with_qual_at(read, pos, qual):
            quals = list(read.query_qualities)
            quals[pos] = qual
            read.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(q + 33) for q in quals)
            )
            return read

        # R1 calls T at pos 1 (an edit), R2 calls the reference C there.
        edited, plain = "ATGTCGATCG", REF_SEQ

        confident_edit = self._pair(edited, plain)
        with_qual_at(confident_edit[0], 1, 40)
        with_qual_at(confident_edit[1], 1, 2)
        assert (
            self._run(tmp_path, fasta_file, confident_edit, min_mapq=0)["editing"][
                "total_edits"
            ]
            == 1
        )

        confident_ref = self._pair(edited, plain)
        with_qual_at(confident_ref[0], 1, 2)
        with_qual_at(confident_ref[1], 1, 40)
        assert (
            self._run(tmp_path, fasta_file, confident_ref, name="b", min_mapq=0)[
                "editing"
            ]["total_edits"]
            == 0
        )

    def test_orphaned_mate_is_still_counted(self, tmp_path, fasta_file):
        # R2 fails --min_mapq, so R1's partner never arrives. It must still be
        # folded in rather than sitting in the buffer and being dropped.
        reads = self._pair(self.EDITED, self.EDITED)
        reads[1].mapping_quality = 0
        m = self._run(tmp_path, fasta_file, reads, min_mapq=20)
        assert m["reads"]["passing"] == 1
        assert m["fragments"]["total"] == 1
        assert m["fragments"]["from_mate_pairs"] == 0
        assert m["fragments"]["from_single_records"] == 1
        assert m["editing"]["total_edits"] == 2

    def test_single_end_reads_are_one_fragment_each(self, tmp_path, fasta_file):
        reads = [_make_read(f"r{i}", self.EDITED, 0, is_paired=False) for i in range(3)]
        m = self._run(tmp_path, fasta_file, reads, min_mapq=0)
        assert m["reads"]["total"] == 3
        assert m["fragments"]["total"] == 3
        assert m["fragments"]["from_mate_pairs"] == 0
        assert m["editing"]["total_opportunities"] == 15  # 5 per fragment

    def test_subsampling_keeps_both_mates_or_neither(self, tmp_path, fasta_file):
        """Sampling is keyed on the read name, so a fragment survives whole.

        Drawing records independently would leave most kept fragments with one
        mate and halve every per-fragment edit count.
        """
        reads = []
        for i in range(500):
            reads += [
                _make_read(
                    f"frag{i}",
                    self.EDITED,
                    0,
                    is_read1=True,
                    mate_pos=0,
                    mate_reverse=True,
                ),
                _make_read(
                    f"frag{i}",
                    self.EDITED,
                    0,
                    is_read1=False,
                    is_reverse=True,
                    mate_pos=0,
                ),
            ]
        m = self._run(tmp_path, fasta_file, reads, min_mapq=0, n_reads=500)

        assert m["sampling"]["subsampled"] is True
        assert 0 < m["fragments"]["total"] < 500
        assert m["fragments"]["from_single_records"] == 0
        assert m["fragments"]["from_mate_pairs"] == m["fragments"]["total"]
        assert m["reads"]["total"] == 2 * m["fragments"]["total"]
        # And the merged rate is the same as in a full run, not half of it.
        assert m["edit_rate_per_fragment"]["mean"] == pytest.approx(2 / 6, abs=1e-6)


class TestTssEnrichment:
    """TSS enrichment as the ENCODE ATAC-seq pipeline defines it.

    The profile is built from insertion 5' ends binned at 10 bp, minus-strand
    windows are flipped, the flanks are normalised to 1, and the score is the
    *peak* of that profile -- not an average over the centre.
    """

    CHROM = "chr1"
    CHROM_LEN = 20_000
    TSS = 10_000
    FLANK = 500  # 100 bins of 10 bp; the outer 10 bins are the background
    READ_LEN = 50

    def _write_bam(self, tmp_path, sites, name="tss.bam"):
        """BAM whose reads have their insertion 5' end at each given position.

        ``sites`` holds ``(position, is_reverse)`` pairs; a reverse read is
        placed so that ``reference_end - 1`` lands on the position.
        """
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"LN": self.CHROM_LEN, "SN": self.CHROM}],
        }
        path = str(tmp_path / name)
        tmp = path + ".unsorted.bam"
        with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
            for i, (pos, is_reverse) in enumerate(sites):
                start = pos - (self.READ_LEN - 1) if is_reverse else pos
                read = _make_read(
                    f"r{i}",
                    "A" * self.READ_LEN,
                    start,
                    is_reverse=is_reverse,
                    is_paired=False,
                )
                bam.write(read)
        pysam.sort("-o", path, tmp)
        os.remove(tmp)
        pysam.index(path)
        return path

    def _bed(self, tmp_path, strand=None, name="tss.bed"):
        path = str(tmp_path / name)
        cols = [self.CHROM, str(self.TSS), str(self.TSS + 1)]
        if strand is not None:
            cols += ["tss1", "0", strand]
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
        return path

    def _flat(self, per_position=1):
        """One insertion at every position of the window: a featureless profile."""
        return [
            (p, False)
            for p in range(self.TSS - self.FLANK, self.TSS + self.FLANK)
            for _ in range(per_position)
        ]

    def _run(self, tmp_path, sites, strand=None):
        bam = self._write_bam(tmp_path, sites)
        bed = self._bed(tmp_path, strand=strand)
        return qc._tss_enrichment(bam, bed, {self.CHROM: self.CHROM_LEN}, 0, self.FLANK)

    def test_flat_signal_scores_one(self, tmp_path):
        res = self._run(tmp_path, self._flat())
        assert res is not None
        assert res.bin_size == 10
        assert len(res.profile) == 2 * self.FLANK // 10
        assert res.n_tss == 1
        # Every bin equals the background, so the peak is the background.
        assert res.score == pytest.approx(1.0)
        assert res.profile.min() == pytest.approx(1.0)

    def test_score_is_the_peak_height(self, tmp_path):
        # Flat background of 1 insertion/bp (10 per bin) plus 90 extra at the
        # TSS, so the centre bin holds 100 against a background of 10.
        sites = self._flat() + [(self.TSS, False)] * 90
        res = self._run(tmp_path, sites)
        assert res is not None
        assert res.background == pytest.approx(10.0)
        assert res.score == pytest.approx(10.0)
        assert res.positions[int(res.profile.argmax())] == pytest.approx(5.0)

    def test_peak_is_found_away_from_the_centre(self, tmp_path):
        # ENCODE takes the maximum of the whole profile. An off-centre peak
        # therefore sets the score; averaging over the centre would miss it.
        sites = self._flat() + [(self.TSS + 200, False)] * 90
        res = self._run(tmp_path, sites)
        assert res is not None
        assert res.score == pytest.approx(10.0)
        assert res.positions[int(res.profile.argmax())] == pytest.approx(205.0)

    def test_reverse_reads_count_their_own_5_prime_end(self, tmp_path):
        # A reverse read's insertion site is reference_end - 1, so these land in
        # the same bin as forward reads starting there.
        sites = self._flat() + [(self.TSS + 200, True)] * 90
        res = self._run(tmp_path, sites)
        assert res is not None
        assert res.positions[int(res.profile.argmax())] == pytest.approx(205.0)

    def test_minus_strand_tss_is_flipped(self, tmp_path):
        # Signal 200 bp to the left of the TSS is *downstream* for a minus-strand
        # gene, so flipping must move the peak to positive coordinates.
        sites = self._flat() + [(self.TSS - 200, False)] * 90
        plus = self._run(tmp_path, sites, strand="+")
        minus = self._run(tmp_path, sites, strand="-")
        assert plus is not None and minus is not None
        assert plus.positions[int(plus.profile.argmax())] < 0
        assert minus.positions[int(minus.profile.argmax())] > 0
        assert minus.score == pytest.approx(plus.score)

    def test_no_usable_tss_returns_none(self, tmp_path):
        bam = self._write_bam(tmp_path, self._flat())
        bed = str(tmp_path / "off.bed")
        with open(bed, "w") as f:
            f.write(f"{self.CHROM}\t10\t11\n")  # window runs off the contig start
        assert (
            qc._tss_enrichment(bam, bed, {self.CHROM: self.CHROM_LEN}, 0, self.FLANK)
            is None
        )

    def test_zero_background_gives_nan_score(self, tmp_path):
        # Insertions only at the TSS: nothing in the flanks to normalise against.
        res = self._run(tmp_path, [(self.TSS, False)] * 20)
        assert res is not None
        assert res.score != res.score  # NaN
        assert res.counts.sum() == 20  # the raw profile is still reported

    def test_flank_smaller_than_one_bin_falls_back_to_1bp(self, tmp_path):
        res = self._run_with_flank(tmp_path, 4)
        assert res is not None
        assert res.bin_size == 1
        assert len(res.profile) == 8

    def _run_with_flank(self, tmp_path, flank):
        bam = self._write_bam(tmp_path, self._flat())
        bed = self._bed(tmp_path)
        return qc._tss_enrichment(bam, bed, {self.CHROM: self.CHROM_LEN}, 0, flank)

    def test_run_qc_writes_csv_and_plots_the_score(self, tmp_path):
        fasta = str(tmp_path / "big.fa")
        with open(fasta, "w") as f:
            f.write(f">{self.CHROM}\n")
            for _ in range(0, self.CHROM_LEN, 60):
                f.write("ACGTCG" * 10 + "\n")
        pysam.faidx(fasta)

        sites = self._flat() + [(self.TSS, False)] * 90
        bam = self._write_bam(tmp_path, sites)
        bed = self._bed(tmp_path)
        out_dir = str(tmp_path / "out")
        m = run_qc(
            bam,
            fasta,
            out_dir,
            "sample",
            tss_path=bed,
            min_mapq=0,
            min_baseq=0,
            tss_flank=self.FLANK,
        )

        tss = m["tss_enrichment"]
        assert tss["score"] == pytest.approx(10.0)
        assert tss["n_tss"] == 1
        assert tss["flank"] == self.FLANK
        assert tss["bin_size"] == 10
        assert tss["profile_csv"] == "sample.tss_enrichment.csv"
        # The vector lives only in the CSV, so the JSON stays small.
        assert "profile" not in tss

        csv_path = os.path.join(out_dir, tss["profile_csv"])
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 2 * self.FLANK // 10
        centre = next(r for r in rows if float(r["position"]) == 5.0)
        assert float(centre["normalized"]) == pytest.approx(10.0, rel=1e-4)
        assert max(float(r["normalized"]) for r in rows) == pytest.approx(tss["score"])
        assert sum(int(r["insertions"]) for r in rows) == tss["total_insertions"]

        html = open(os.path.join(out_dir, "sample.html")).read()
        assert "TSS enrichment" in html
        assert "sample.tss_enrichment.csv" in html
        # Its own plot, on top of the multi-panel summary figure.
        assert html.count("data:image/png;base64,") >= 3


class TestSubsampling:
    """`n_reads` draws a uniform sample instead of reading everything."""

    def _bam_with(self, tmp_path, n_reads, fasta_file):
        """A BAM of `n_reads` identical reads, each carrying one C->T edit.

        Identical reads make the editing rate independent of which subset is
        drawn, so a rate that changes under sampling is a real bug rather than
        sampling noise.
        """
        # REF_SEQ is ACGTCGATCG; flip the C at index 1 to T.
        edited = REF_SEQ[:1] + "T" + REF_SEQ[2:]
        reads = [_make_read(f"r{i}", edited, 0) for i in range(n_reads)]
        return _write_bam(str(tmp_path / "s.bam"), reads)

    def test_none_uses_every_read(self, tmp_path, fasta_file):
        bam = self._bam_with(tmp_path, 200, fasta_file)
        m = run_qc(bam, fasta_file, str(tmp_path / "o"), "all", plot=False)
        assert m["reads"]["total"] == 200
        assert m["sampling"]["subsampled"] is False
        assert m["sampling"]["fraction"] == 1.0

    def test_subsamples_to_about_the_requested_count(self, tmp_path, fasta_file):
        bam = self._bam_with(tmp_path, 400, fasta_file)
        m = run_qc(bam, fasta_file, str(tmp_path / "o"), "sub", plot=False, n_reads=100)
        assert m["sampling"]["subsampled"] is True
        assert m["sampling"]["fraction"] == pytest.approx(0.25)
        # Binomial(400, 0.25): sd = 8.7, so a wide band is still a real check.
        assert 50 <= m["reads"]["total"] <= 160
        assert m["reads"]["total"] < 400

    def test_is_reproducible(self, tmp_path, fasta_file):
        bam = self._bam_with(tmp_path, 300, fasta_file)
        a = run_qc(bam, fasta_file, str(tmp_path / "a"), "x", plot=False, n_reads=90)
        b = run_qc(bam, fasta_file, str(tmp_path / "b"), "x", plot=False, n_reads=90)
        assert a["reads"]["total"] == b["reads"]["total"]
        assert a["editing"] == b["editing"]

    def test_larger_request_than_the_bam_uses_all_reads(self, tmp_path, fasta_file):
        bam = self._bam_with(tmp_path, 50, fasta_file)
        m = run_qc(
            bam, fasta_file, str(tmp_path / "o"), "big", plot=False, n_reads=10_000
        )
        assert m["reads"]["total"] == 50
        assert m["sampling"]["subsampled"] is False

    def test_rejects_non_positive(self, tmp_path, fasta_file):
        bam = self._bam_with(tmp_path, 20, fasta_file)
        with pytest.raises(ValueError, match="n_reads must be positive"):
            run_qc(bam, fasta_file, str(tmp_path / "o"), "z", plot=False, n_reads=0)

    def test_draw_ignores_whether_a_read_carries_an_edit(self, tmp_path, fasta_file):
        """The coin is flipped on the read, never on its content.

        Every other read here carries no editing event at all. A sampler that
        required an edit -- or merely favoured edited reads -- would push
        mean_edits_per_read up towards 1; a blind draw leaves it at 0.5.
        """
        edited = REF_SEQ[:1] + "T" + REF_SEQ[2:]  # one C->T at position 1
        reads = [
            _make_read(f"r{i}", edited if i % 2 else REF_SEQ, 0) for i in range(1000)
        ]
        bam = _write_bam(str(tmp_path / "mix.bam"), reads)

        full = run_qc(bam, fasta_file, str(tmp_path / "f"), "f", plot=False)
        sub = run_qc(bam, fasta_file, str(tmp_path / "s"), "s", plot=False, n_reads=500)

        assert full["editing"]["mean_edits_per_fragment"] == pytest.approx(0.5)
        # Binomial(~500, 0.5) on the sample: sd ~= 0.022, so this band is ~5 sd.
        assert sub["editing"]["mean_edits_per_fragment"] == pytest.approx(0.5, abs=0.11)
        assert 0 < sub["reads"]["total"] < 1000

    def test_editing_rate_is_unbiased_by_sampling(self, tmp_path, fasta_file):
        """A rate must survive subsampling; only absolute counts shrink."""
        bam = self._bam_with(tmp_path, 600, fasta_file)
        full = run_qc(bam, fasta_file, str(tmp_path / "f"), "f", plot=False)
        sub = run_qc(bam, fasta_file, str(tmp_path / "s"), "s", plot=False, n_reads=300)
        assert sub["reads"]["total"] < full["reads"]["total"]
        assert sub["editing"]["global_edit_rate"] == pytest.approx(
            full["editing"]["global_edit_rate"], abs=1e-9
        )
        assert sub["editing"]["global_edit_rate"] > 0  # the rate is real, not 0 == 0
