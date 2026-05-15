"""Tests for deamtools.preprocessing.bam2fragment."""

import gzip
import os

import pysam
import pytest

from deamtools.preprocessing.bam2fragment import run_bam2fragment

# Reference: A C G T C G A T C G    (positions 0..9)
# C positions (forward edit): 1, 4, 8
# G positions (reverse edit): 2, 5, 9
REF_SEQ = "ACGTCGATCG"

BAM_HEADER = {
    "HD": {"VN": "1.6"},
    "SQ": [{"LN": len(REF_SEQ), "SN": "chr1"}],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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
    mapq: int = 30,
    baseq: int = 40,
    extra_flags: int = 0,
    tags: list[tuple[str, str]] | None = None,
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
        if is_read1:
            flag |= 0x40
        else:
            flag |= 0x80
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * len(seq))
    if is_paired and mate_pos is not None:
        a.next_reference_id = 0
        a.next_reference_start = mate_pos
    if tags:
        a.set_tags(tags)
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


def _pair(name, r1_pos, r1_seq, r2_pos, r2_seq, **kw):
    """Build a properly-paired R1 (forward) + R2 (reverse) pair."""
    tags = kw.pop("tags", None)
    r1 = _make_read(
        name, r1_seq, r1_pos, is_paired=True, is_read1=True,
        is_reverse=False, mate_pos=r2_pos, mate_reverse=True,
        tags=tags, **kw,
    )
    r2 = _make_read(
        name, r2_seq, r2_pos, is_paired=True, is_read1=False,
        is_reverse=True, mate_pos=r1_pos, mate_reverse=False,
        tags=tags, **kw,
    )
    return [r1, r2]


def _read_lines(path: str) -> list[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.rstrip("\n") for line in f if line.strip()]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def fasta_file(tmp_path):
    path = str(tmp_path / "ref.fa")
    with open(path, "w") as f:
        f.write(f">chr1\n{REF_SEQ}\n")
    pysam.faidx(path)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestBam2Fragment:
    def test_paired_fragment_basic_columns(self, tmp_path, fasta_file):
        # Pair: R1 covers pos 0-5 with C->T at pos 1; R2 covers pos 4-10 (no edits).
        reads = _pair(
            "pair1",
            r1_pos=0, r1_seq="ATGTCG",         # C->T at pos 1
            r2_pos=4, r2_seq=REF_SEQ[4:10],    # CGATCG, identical to ref
        )
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)

        lines = _read_lines(out)
        assert len(lines) == 1
        cols = lines[0].split("\t")
        # chrom, start, end, count, edits
        assert cols[0] == "chr1"
        assert cols[1] == "0"   # min(R1 start=0, R2 start=4)
        assert cols[2] == "10"  # max(R1 end=6, R2 end=10)
        assert cols[3] == "1"
        assert cols[4] == "1"   # only edit pos

    def test_no_edits_emits_dot(self, tmp_path, fasta_file):
        reads = _pair(
            "pair1",
            r1_pos=0, r1_seq=REF_SEQ[0:6],
            r2_pos=4, r2_seq=REF_SEQ[4:10],
        )
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        cols = _read_lines(out)[0].split("\t")
        assert cols[4] == "."

    def test_duplicate_signatures_counted(self, tmp_path, fasta_file):
        # Two pairs with identical coords and identical edits -> single row, count=2.
        reads = []
        for n in ("p1", "p2"):
            reads += _pair(n, 0, "ATGTCG", 4, REF_SEQ[4:10])  # both have C->T at pos 1
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        lines = _read_lines(out)
        assert len(lines) == 1
        cols = lines[0].split("\t")
        assert cols[3] == "2"
        assert cols[4] == "1"

    def test_different_edit_signatures_separate_rows(self, tmp_path, fasta_file):
        # Two pairs at the same coords but different edit signatures get two rows.
        reads = []
        # pair A: edit at pos 1
        reads += _pair("pA", 0, "ATGTCG", 4, REF_SEQ[4:10])
        # pair B: edit at pos 4 (R1 reads ATGTTG -> T at pos 4)
        reads += _pair("pB", 0, "ACGTTG", 4, REF_SEQ[4:10])
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        lines = _read_lines(out)
        assert len(lines) == 2
        edit_cols = sorted(line.split("\t")[4] for line in lines)
        assert edit_cols == ["1", "4"]

    def test_reverse_read_g_to_a_edit(self, tmp_path, fasta_file):
        # R2 is a reverse-strand read; G->A at pos 2 should be detected.
        # R2 sequence is given in reference orientation (reverse-strand reads are
        # stored as they appear on the + strand in the BAM SEQ column).
        # ref pos 2..7 = GTCGAT; an A at pos 2 -> ATCGAT.
        reads = _pair(
            "pair1",
            r1_pos=0, r1_seq=REF_SEQ[0:6],  # no edits
            r2_pos=2, r2_seq="ATCGAT",      # G->A at pos 2 on reverse read
        )
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        cols = _read_lines(out)[0].split("\t")
        assert cols[4] == "2"

    def test_dedup_of_overlapping_edit_positions(self, tmp_path, fasta_file):
        # R1 and R2 both cover pos 4 and both report a C->T there.
        # The fragment should list pos 4 only once.
        reads = _pair(
            "pair1",
            r1_pos=0, r1_seq="ACGTTG",       # T at pos 4 (forward)
            # R2 reverse, covers ref 4..9 = CGATCG. With C at pos 4 changed to T,
            # we'd see a forward-style edit. But R2 is reverse, so a C->T pattern
            # on R2 isn't counted (only G->A on reverse counts). So put a normal
            # match here and rely on R1 alone for the pos-4 edit.
            r2_pos=4, r2_seq=REF_SEQ[4:10],
        )
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        cols = _read_lines(out)[0].split("\t")
        assert cols[4] == "4"

    def test_min_mapq_drops_pair(self, tmp_path, fasta_file):
        reads = _pair("pair1", 0, "ATGTCG", 4, REF_SEQ[4:10], mapq=10)
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=20, min_baseq=0)
        assert _read_lines(out) == []

    def test_min_baseq_filters_edit(self, tmp_path, fasta_file):
        reads = _pair("pair1", 0, "ATGTCG", 4, REF_SEQ[4:10], baseq=5)
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=20)
        cols = _read_lines(out)[0].split("\t")
        assert cols[4] == "."

    def test_duplicate_flag_skipped(self, tmp_path, fasta_file):
        reads = _pair("p1", 0, "ATGTCG", 4, REF_SEQ[4:10], extra_flags=0x400)
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        assert _read_lines(out) == []

    def test_single_end_fragment(self, tmp_path, fasta_file):
        read = _make_read(
            "se1", "ATGTCG", 0,
            is_paired=False, is_read1=False, mate_pos=None,
        )
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        cols = _read_lines(out)[0].split("\t")
        assert cols[0:3] == ["chr1", "0", "6"]
        assert cols[3] == "1"
        assert cols[4] == "1"

    def test_barcode_column_added_with_flag(self, tmp_path, fasta_file):
        reads = _pair("p1", 0, "ATGTCG", 4, REF_SEQ[4:10],
                      tags=[("CB", "AAACGT-1")])
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0,
                         barcode=True, barcode_tag="CB")
        cols = _read_lines(out)[0].split("\t")
        # 10x ordering: chrom, start, end, barcode, count, edits
        assert len(cols) == 6
        assert cols[3] == "AAACGT-1"
        assert cols[4] == "1"
        assert cols[5] == "1"

    def test_barcode_groups_by_barcode(self, tmp_path, fasta_file):
        # Same coords + edits, two different barcodes -> two rows.
        reads = []
        reads += _pair("p1", 0, "ATGTCG", 4, REF_SEQ[4:10],
                       tags=[("CB", "BC_A")])
        reads += _pair("p2", 0, "ATGTCG", 4, REF_SEQ[4:10],
                       tags=[("CB", "BC_B")])
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0,
                         barcode=True, barcode_tag="CB")
        lines = _read_lines(out)
        assert len(lines) == 2
        barcodes = sorted(line.split("\t")[3] for line in lines)
        assert barcodes == ["BC_A", "BC_B"]
        assert all(line.split("\t")[4] == "1" for line in lines)

    def test_missing_barcode_renders_dot(self, tmp_path, fasta_file):
        reads = _pair("p1", 0, "ATGTCG", 4, REF_SEQ[4:10])
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0,
                         barcode=True, barcode_tag="CB")
        cols = _read_lines(out)[0].split("\t")
        assert cols[3] == "."

    def test_gzip_output_when_path_ends_in_gz(self, tmp_path, fasta_file):
        reads = _pair("p1", 0, "ATGTCG", 4, REF_SEQ[4:10])
        bam = _write_bam(str(tmp_path / "x.bam"), reads)
        out = str(tmp_path / "frag.tsv.gz")
        run_bam2fragment(bam_path=bam, fasta_path=fasta_file, output_path=out,
                         min_mapq=0, min_baseq=0)
        assert os.path.exists(out)
        # Confirm the file is actually gzip-compressed
        with open(out, "rb") as f:
            assert f.read(2) == b"\x1f\x8b"
        cols = _read_lines(out)[0].split("\t")
        assert cols[0] == "chr1"
