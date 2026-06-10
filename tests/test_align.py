"""Tests for deamtools.align.align (pure-Python conversion / restoration).

These do not require ``bwa`` or ``samtools``; they exercise the read converter
and the SAM post-processing that restores original sequences and chromosome
names.
"""

import io
from types import SimpleNamespace

import pytest

from deamtools.align.align import (
    CT_TABLE,
    GA_TABLE,
    _emit_clean_header,
    _feed_converted,
    _hard_clip_offsets,
    _process_sam,
    _restore_alignment,
    _revcomp,
    _write_converted,
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


class TestWriteConverted:
    def test_ct_conversion_and_ys_tag(self):
        out = io.StringIO()
        entry = SimpleNamespace(name="r1", sequence="ACGTC", quality="IIIII")
        _write_converted(out, entry, CT_TABLE)
        assert out.getvalue() == "@r1\tYS:Z:ACGTC\nATGTT\n+\nIIIII\n"

    def test_ga_conversion(self):
        out = io.StringIO()
        entry = SimpleNamespace(name="r2", sequence="ACGTC", quality="IIIII")
        _write_converted(out, entry, GA_TABLE)
        assert out.getvalue() == "@r2\tYS:Z:ACGTC\nACATC\n+\nIIIII\n"

    def test_missing_quality_filled(self):
        out = io.StringIO()
        entry = SimpleNamespace(name="r1", sequence="ACGTC", quality=None)
        _write_converted(out, entry, CT_TABLE)
        assert out.getvalue() == "@r1\tYS:Z:ACGTC\nATGTT\n+\nIIIII\n"


class TestFeedConverted:
    def test_single_end_converts_ct(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        _write_fastq(fq1, [("a", "ACGTC", "IIIII"), ("b", "CCGG", "IIII")])
        cap = _Capture()
        _feed_converted(fq1, None, cap)
        text = cap.getvalue()
        # Both records present, C->T converted.
        assert "@a\tYS:Z:ACGTC\nATGTT\n" in text
        assert "@b\tYS:Z:CCGG\nTTGG\n" in text

    def test_paired_interleaves_r1_ct_r2_ga(self, tmp_path):
        fq1 = str(tmp_path / "r1.fq")
        fq2 = str(tmp_path / "r2.fq")
        _write_fastq(fq1, [("a", "ACGTC", "IIIII")])
        _write_fastq(fq2, [("a", "ACGTC", "IIIII")])
        cap = _Capture()
        _feed_converted(fq1, fq2, cap)
        # R1 is C->T, R2 is G->A, in that order.
        assert cap.getvalue() == (
            "@a\tYS:Z:ACGTC\nATGTT\n+\nIIIII\n"
            "@a\tYS:Z:ACGTC\nACATC\n+\nIIIII\n"
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
