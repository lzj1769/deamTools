"""Tests for deamtools.align.align (pure-Python conversion / restoration).

These do not require ``bwa`` or ``samtools``; they exercise the read converter
and the SAM post-processing that restores original sequences and chromosome
names.
"""

import io

import pytest

from deamtools.align.align import (
    _emit_clean_header,
    _feed_converted,
    _hard_clip_offsets,
    _process_sam,
    _restore_alignment,
    _revcomp,
    _write_record,
    run_align,
)


class _Capture(io.StringIO):
    """StringIO whose close() is a no-op so content survives _feed_converted."""

    def close(self) -> None:  # noqa: D401
        pass


def _write_fastq(path: str, reads: list[tuple[str, str, str]]) -> None:
    with open(path, "w") as f:
        for name, seq, qual in reads:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


class TestHelpers:
    def test_revcomp(self):
        assert _revcomp("ACGTN") == "NACGT"
        assert _revcomp("acgt") == "acgt"  # palindrome, case preserved

    def test_hard_clip_offsets(self):
        assert _hard_clip_offsets("10H50M5H") == (10, 5)
        assert _hard_clip_offsets("50M") == (0, 0)
        assert _hard_clip_offsets("5S50M") == (0, 0)  # soft clips are not hard
        assert _hard_clip_offsets("3H50M") == (3, 0)


class TestWriteRecord:
    def test_record_with_ys_and_yc_tags(self):
        out = io.StringIO()
        _write_record(out, "r1", "ATGTT", "IIIII", "ct", "ACGTC")
        assert out.getvalue() == "@r1\tYS:Z:ACGTC\tYC:Z:ct\nATGTT\n+\nIIIII\n"


class TestFeedConverted:
    def test_single_end_emits_ct_and_ga_candidates(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        _write_fastq(fq1, [("a", "ACGTC", "IIIII")])
        cap = _Capture()
        _feed_converted(fq1, None, cap)
        # Two candidates per read: C->T (ct) then G->A (ga), same read name.
        assert cap.getvalue() == (
            "@a\tYS:Z:ACGTC\tYC:Z:ct\nATGTT\n+\nIIIII\n"
            "@a\tYS:Z:ACGTC\tYC:Z:ga\nACATC\n+\nIIIII\n"
        )

    def test_missing_quality_filled(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        with open(fq1, "w") as f:  # FASTA-style (no qualities) via FastxFile
            f.write(">a\nACGTC\n")
        cap = _Capture()
        _feed_converted(fq1, None, cap)
        assert "\nATGTT\n+\nIIIII\n" in cap.getvalue()

    def test_paired_emits_f_and_r_orientations(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        fq2 = str(tmp_path / "r2.fq")
        _write_fastq(fq1, [("a", "ACGTC", "IIIII")])
        _write_fastq(fq2, [("a", "ACGTC", "JJJJJ")])
        cap = _Capture()
        _feed_converted(fq1, fq2, cap)
        assert cap.getvalue() == (
            # orientation f: r1 C->T, r2 G->A
            "@a\tYS:Z:ACGTC\tYC:Z:f\nATGTT\n+\nIIIII\n"
            "@a\tYS:Z:ACGTC\tYC:Z:f\nACATC\n+\nJJJJJ\n"
            # orientation r: r1 G->A, r2 C->T
            "@a\tYS:Z:ACGTC\tYC:Z:r\nACATC\n+\nIIIII\n"
            "@a\tYS:Z:ACGTC\tYC:Z:r\nATGTT\n+\nJJJJJ\n"
        )

    def test_paired_length_mismatch_raises(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        fq2 = str(tmp_path / "r2.fq")
        _write_fastq(fq1, [("a", "ACGT", "IIII"), ("b", "ACGT", "IIII")])
        _write_fastq(fq2, [("a", "ACGT", "IIII")])
        with pytest.raises(ValueError):
            _feed_converted(fq1, fq2, _Capture())


def _sam_line(fields: list[str]) -> str:
    return "\t".join(fields) + "\n"


class TestRestoreAlignment:
    def test_strips_f_r_prefix_from_rname_and_rnext(self):
        line = _sam_line(
            ["r1", "0", "fchr1", "10", "60", "5M", "rchr1", "20", "0",
             "ATGTT", "IIIII", "YS:Z:ACGTC"]
        )
        out = _restore_alignment(line).rstrip("\n").split("\t")
        assert out[2] == "chr1"  # RNAME prefix stripped
        assert out[6] == "chr1"  # RNEXT prefix stripped

    def test_restores_original_seq_forward(self):
        line = _sam_line(
            ["r1", "0", "fchr1", "10", "60", "5M", "*", "0", "0",
             "ATGTT", "IIIII", "YS:Z:ACGTC", "NM:i:1"]
        )
        out = _restore_alignment(line)
        fields = out.rstrip("\n").split("\t")
        assert fields[9] == "ACGTC"      # SEQ restored from YS
        assert "YS:Z:" not in out        # YS tag dropped
        assert "NM:i:1" in fields        # other tags kept

    def test_restores_revcomp_seq_on_reverse_strand(self):
        # Reverse-strand read: BAM SEQ is revcomp of the converted read; the
        # restored SEQ must be revcomp of the original (YS) read.
        line = _sam_line(
            ["r1", "16", "fchr1", "10", "60", "5M", "*", "0", "0",
             "AACAT", "IIIII", "YS:Z:ACGTC"]
        )
        fields = _restore_alignment(line).rstrip("\n").split("\t")
        assert fields[9] == _revcomp("ACGTC")  # "GACGT"

    def test_hard_clip_trims_original_seq(self):
        # CIGAR 2H5M: SEQ is hard-clipped to 5 bp; restored SEQ is YS[2:7].
        line = _sam_line(
            ["r1", "0", "chr1", "10", "60", "2H5M", "*", "0", "0",
             "XXXXX", "IIIII", "YS:Z:ACGTCGA"]
        )
        fields = _restore_alignment(line).rstrip("\n").split("\t")
        assert fields[9] == "GTCGA"  # ACGTCGA[2:7]

    def test_short_line_passthrough(self):
        line = "not\ta\tsam\trecord\n"
        assert _restore_alignment(line) == line

    def test_unmapped_rname_star_not_stripped(self):
        line = _sam_line(
            ["r1", "4", "*", "0", "0", "*", "*", "0", "0",
             "ATGTT", "IIIII", "YS:Z:ACGTC"]
        )
        fields = _restore_alignment(line).rstrip("\n").split("\t")
        assert fields[2] == "*"
        assert fields[9] == "ACGTC"  # SEQ still restored

    def test_drops_yc_tag(self):
        line = _sam_line(
            ["r1", "0", "fchr1", "10", "60", "5M", "*", "0", "0",
             "ATGTT", "IIIII", "YS:Z:ACGTC", "YC:Z:f", "NM:i:1"]
        )
        out = _restore_alignment(line)
        assert "YC:Z:" not in out          # candidate marker dropped
        assert "YS:Z:" not in out          # original-seq tag dropped
        assert "NM:i:1" in out             # real tags kept
        assert out.rstrip("\n").split("\t")[9] == "ACGTC"


class TestIndexResolution:
    """run_align resolves the converted-index location before touching bwa.

    The missing-``.bwt`` check happens first, so these exercise index-path
    resolution without needing bwa/samtools installed.
    """

    def test_default_index_is_next_to_fasta(self, tmp_path):
        fasta = str(tmp_path / "ref.fa")
        with open(fasta, "w") as f:
            f.write(">c\nACGT\n")
        with pytest.raises(FileNotFoundError, match=r"ref\.fa\.deamtools\.c2t\.bwt"):
            run_align(fasta, "r1.fq", str(tmp_path), "out")

    def test_custom_index_path_is_used(self, tmp_path):
        fasta = str(tmp_path / "ref.fa")
        with open(fasta, "w") as f:
            f.write(">c\nACGT\n")
        custom = str(tmp_path / "idx" / "myref.deamtools.c2t")
        with pytest.raises(FileNotFoundError, match=r"myref\.deamtools\.c2t\.bwt"):
            run_align(
                fasta, "r1.fq", str(tmp_path), "out", index_path=custom
            )

    def test_custom_index_path_found_proceeds_past_bwt_check(self, tmp_path):
        # With a present .bwt at the custom location, the .bwt check passes and
        # run_align advances to the FASTA-index check (proving --index is used).
        fasta = str(tmp_path / "ref.fa")
        with open(fasta, "w") as f:
            f.write(">c\nACGT\n")
        idx_dir = tmp_path / "idx"
        idx_dir.mkdir()
        custom = str(idx_dir / "myref.deamtools.c2t")
        open(custom + ".bwt", "w").close()  # satisfy the index check
        # No .fai next to the FASTA -> next check fails, mentioning .fai.
        with pytest.raises(FileNotFoundError, match=r"ref\.fa\.fai"):
            run_align(fasta, "r1.fq", str(tmp_path), "out", index_path=custom)


class TestHeaderAndStream:
    def test_emit_clean_header(self, tmp_path):
        fasta = str(tmp_path / "ref.fa")
        with open(fasta + ".fai", "w") as f:
            f.write("chr1\t1000\t6\t60\t61\n")
            f.write("chr2\t2000\t1100\t60\t61\n")
        out = io.StringIO()
        _emit_clean_header(fasta, out)
        assert out.getvalue() == (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "@SQ\tSN:chr2\tLN:2000\n"
        )

    def test_process_sam_drops_hd_sq_keeps_pg_and_restores(self):
        bwa_out = io.StringIO(
            "@HD\tVN:1.6\n"
            "@SQ\tSN:fchr1\tLN:1000\n"
            "@PG\tID:bwa\tPN:bwa\n"
            + _sam_line(
                ["r1", "0", "fchr1", "10", "60", "5M", "*", "0", "0",
                 "ATGTT", "IIIII", "YS:Z:ACGTC"]
            )
        )
        sink = io.StringIO()
        _process_sam(bwa_out, sink)
        text = sink.getvalue()
        assert "@HD" not in text          # dropped (re-emitted elsewhere)
        assert "@SQ" not in text          # dropped (f/r contigs)
        assert "@PG\tID:bwa" in text      # kept
        # Alignment line restored: prefix stripped + SEQ recovered.
        aln = text.strip().split("\n")[-1].split("\t")
        assert aln[2] == "chr1"
        assert aln[9] == "ACGTC"

    def test_process_sam_se_picks_higher_scoring_candidate(self):
        # Same read name, two candidates; ct has higher AS -> ct kept, ga dropped.
        ct = _sam_line(
            ["r1", "0", "fchr1", "10", "60", "5M", "*", "0", "0",
             "ATGTT", "IIIII", "YS:Z:ACGTC", "YC:Z:ct", "AS:i:50"]
        )
        ga = _sam_line(
            ["r1", "0", "rchr2", "20", "60", "5M", "*", "0", "0",
             "ACATC", "IIIII", "YS:Z:ACGTC", "YC:Z:ga", "AS:i:10"]
        )
        sink = io.StringIO()
        _process_sam(io.StringIO(ct + ga), sink)
        lines = [ln for ln in sink.getvalue().splitlines() if ln]
        assert len(lines) == 1                 # only the winning candidate
        fields = lines[0].split("\t")
        assert fields[2] == "chr1"             # ct candidate (fchr1 -> chr1)
        assert "YC:Z:" not in lines[0]
        assert fields[9] == "ACGTC"

    def test_process_sam_pe_picks_best_orientation(self):
        # Orientation r (40+40) beats f (10+10); both mates of r are kept.
        f1 = _sam_line(["p", "65", "fchr1", "10", "60", "5M", "=", "30", "25",
                        "ATGTT", "IIIII", "YS:Z:ACGTC", "YC:Z:f", "AS:i:10"])
        f2 = _sam_line(["p", "129", "fchr1", "30", "60", "5M", "=", "10", "-25",
                        "ACATC", "IIIII", "YS:Z:ACGTC", "YC:Z:f", "AS:i:10"])
        r1 = _sam_line(["p", "65", "rchr1", "10", "60", "5M", "=", "30", "25",
                        "ACATC", "IIIII", "YS:Z:ACGTC", "YC:Z:r", "AS:i:40"])
        r2 = _sam_line(["p", "129", "rchr1", "30", "60", "5M", "=", "10", "-25",
                        "ATGTT", "IIIII", "YS:Z:ACGTC", "YC:Z:r", "AS:i:40"])
        sink = io.StringIO()
        _process_sam(io.StringIO(f1 + f2 + r1 + r2), sink)
        lines = [ln for ln in sink.getvalue().splitlines() if ln]
        assert len(lines) == 2                 # both mates of the winner only
        for ln in lines:
            assert ln.split("\t")[2] == "chr1"  # rchr1 -> chr1
            assert "YC:Z:" not in ln
