"""Tests for deamtools.motif.match (motifmatchpy scanning).

These build motifs in memory (no JASPAR/pyjaspar needed); motifmatchpy is a
hard dependency of the package.
"""

import os
from types import SimpleNamespace

import pysam
import pytest

from deamtools.motif.match import (
    load_motifs_from_files,
    prepare_scanner,
    run_motif_matching,
    scan_sequence,
)


def _motif(consensus: str, name: str = "TEST", peak: int = 100):
    """A deterministic motif whose consensus is `consensus`."""
    counts = {b: [0] * len(consensus) for b in "ACGT"}
    for j, base in enumerate(consensus):
        counts[base][j] = peak
    return SimpleNamespace(matrix_id=name, name=name, counts=counts)


class TestScanSequence:
    def test_forward_match_coordinates(self):
        # Consensus AAACCC (not its own reverse complement: rc = GGGTTT).
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        # AAACCC sits at offset 2 within the sequence.
        matches = scan_sequence(scanner, "TTAAACCCTT", "chr1", offset=0)
        plus = [m for m in matches if m[5] == "+"]
        assert len(plus) == 1
        chrom, start, end, name, score, strand = plus[0]
        assert (chrom, start, end, strand) == ("chr1", 2, 8, "+")
        assert name == "TEST"
        assert score > 0

    def test_offset_is_added_to_position(self):
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, "TTAAACCCTT", "chr1", offset=1000)
        plus = [m for m in matches if m[5] == "+"]
        assert plus[0][1] == 1002 and plus[0][2] == 1008

    def test_reverse_strand_match(self):
        # The reverse complement of AAACCC is GGGTTT; placing GGGTTT in the
        # sequence yields a minus-strand hit.
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, "TTGGGTTTTT", "chr1")
        minus = [m for m in matches if m[5] == "-"]
        assert len(minus) == 1
        assert (minus[0][1], minus[0][2]) == (2, 8)

    def test_no_match_returns_empty(self):
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, "TTTTTTTTTT", "chr1")
        assert matches == []


class TestRunMotifMatching:
    def _fasta(self, tmp_path, seq):
        path = str(tmp_path / "ref.fa")
        with open(path, "w") as f:
            f.write(f">chr1\n{seq}\n")
        pysam.faidx(path)
        return path

    def test_writes_bed_of_matches(self, tmp_path):
        # AAACCC at genomic position 5.
        seq = "TTTTT" + "AAACCC" + "TTTTT"  # len 16, motif at [5, 11)
        fasta = self._fasta(tmp_path, seq)
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t0\t16\n")
        out = str(tmp_path / "mpbs.bed")

        run_motif_matching(
            fasta,
            str(bed),
            str(tmp_path),
            "mpbs",
            motifs=[_motif("AAACCC")],
            p_value=1e-3,
        )

        lines = [ln for ln in open(out).read().splitlines() if ln]
        plus = [ln.split("\t") for ln in lines if ln.endswith("+")]
        assert len(plus) == 1
        cols = plus[0]
        assert cols[0] == "chr1"
        assert cols[1] == "5" and cols[2] == "11"
        assert cols[3] == "TEST"
        assert cols[5] == "+"

    def test_offset_within_region(self, tmp_path):
        # Region starts at 100; the motif is 3 bp into the fetched window.
        seq = "A" * 100 + "TTT" + "AAACCC" + "TTT"  # motif at genomic 103
        fasta = self._fasta(tmp_path, seq)
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t100\t115\n")
        out = str(tmp_path / "mpbs.bed")

        run_motif_matching(
            fasta,
            str(bed),
            str(tmp_path),
            "mpbs",
            motifs=[_motif("AAACCC")],
            p_value=1e-3,
        )
        plus = [
            ln.split("\t") for ln in open(out).read().splitlines() if ln.endswith("+")
        ]
        assert plus and plus[0][1] == "103" and plus[0][2] == "109"

    def test_nested_output_dir_created(self, tmp_path):
        fasta = self._fasta(tmp_path, "TTTTTAAACCCTTTTT")
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t0\t16\n")
        out = str(tmp_path / "sub" / "dir" / "mpbs.bed")
        run_motif_matching(
            fasta,
            str(bed),
            str(tmp_path / "sub" / "dir"),
            "mpbs",
            motifs=[_motif("AAACCC")],
            p_value=1e-3,
        )
        assert os.path.exists(out)


def _write_pfm(path, consensus, peak=100):
    """A PFM whose consensus is `consensus`, in the 4-row A/C/G/T layout."""
    rows = []
    for base in "ACGT":
        rows.append(" ".join(str(peak if b == base else 0) for b in consensus))
    path.write_text("\n".join(rows) + "\n")
    return str(path)


class TestMotifFiles:
    def test_name_comes_from_the_file_stem(self, tmp_path):
        path = _write_pfm(tmp_path / "MA0001.1.pfm", "AAACCC")
        (motif,) = load_motifs_from_files([path])
        assert motif.name == "MA0001.1"
        assert motif.length == 6

    def test_missing_file_is_reported_before_parsing(self, tmp_path):
        good = _write_pfm(tmp_path / "ok.pfm", "AAACCC")
        with pytest.raises(FileNotFoundError, match="nope.pfm"):
            load_motifs_from_files([good, str(tmp_path / "nope.pfm")])

    def test_file_motifs_scan_like_count_motifs(self, tmp_path):
        """A PFM file and the equivalent in-memory motif must agree."""
        path = _write_pfm(tmp_path / "AAACCC.pfm", "AAACCC")
        from_file = prepare_scanner(load_motifs_from_files([path]), p_value=1e-3)
        from_counts = prepare_scanner([_motif("AAACCC", name="AAACCC")], p_value=1e-3)
        a = scan_sequence(from_file, "TTAAACCCTT", "chr1")
        b = scan_sequence(from_counts, "TTAAACCCTT", "chr1")
        assert [(x[1], x[2], x[5]) for x in a] == [(x[1], x[2], x[5]) for x in b]

    def test_run_motif_matching_accepts_motif_files(self, tmp_path):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\n" + "TT" + "AAACCC" + "TT" * 20 + "\n")
        pysam.faidx(str(fasta))
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t48\n")
        path = _write_pfm(tmp_path / "EBOX.pfm", "AAACCC")

        run_motif_matching(
            fasta_path=str(fasta),
            bed_path=str(bed),
            out_dir=str(tmp_path / "out"),
            out_name="hits",
            motif_files=[path],
            p_value=1e-3,
        )
        lines = (tmp_path / "out" / "hits.bed").read_text().splitlines()
        assert lines
        assert all(line.split("\t")[3] == "EBOX" for line in lines)

    def test_motif_files_take_precedence_over_motifs(self, tmp_path):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\n" + "TTAAACCCTT" + "\n")
        pysam.faidx(str(fasta))
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t10\n")
        path = _write_pfm(tmp_path / "FROMFILE.pfm", "AAACCC")

        run_motif_matching(
            fasta_path=str(fasta),
            bed_path=str(bed),
            out_dir=str(tmp_path / "out"),
            out_name="hits",
            motifs=[_motif("AAACCC", name="FROMMEMORY")],
            motif_files=[path],
            p_value=1e-3,
        )
        names = {
            line.split("\t")[3]
            for line in (tmp_path / "out" / "hits.bed").read_text().splitlines()
        }
        assert names == {"FROMFILE"}
