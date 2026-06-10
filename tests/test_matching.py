"""Tests for deamtools.motif.matching (MOODS scanning).

These build motifs in memory (no JASPAR/pyjaspar needed); MOODS is a hard
dependency of the package.
"""

import os
from types import SimpleNamespace

import pysam

from deamtools.motif.matching import (
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
        matches = scan_sequence(scanner, [motif], "TTAAACCCTT", "chr1", offset=0)
        plus = [m for m in matches if m[5] == "+"]
        assert len(plus) == 1
        chrom, start, end, name, score, strand = plus[0]
        assert (chrom, start, end, strand) == ("chr1", 2, 8, "+")
        assert name == "TEST"
        assert score > 0

    def test_offset_is_added_to_position(self):
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, [motif], "TTAAACCCTT", "chr1", offset=1000)
        plus = [m for m in matches if m[5] == "+"]
        assert plus[0][1] == 1002 and plus[0][2] == 1008

    def test_reverse_strand_match(self):
        # The reverse complement of AAACCC is GGGTTT; placing GGGTTT in the
        # sequence yields a minus-strand hit.
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, [motif], "TTGGGTTTTT", "chr1")
        minus = [m for m in matches if m[5] == "-"]
        assert len(minus) == 1
        assert (minus[0][1], minus[0][2]) == (2, 8)

    def test_no_match_returns_empty(self):
        motif = _motif("AAACCC")
        scanner = prepare_scanner([motif], p_value=1e-3)
        matches = scan_sequence(scanner, [motif], "TTTTTTTTTT", "chr1")
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
            fasta, str(bed), out, motifs=[_motif("AAACCC")], p_value=1e-3
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
            fasta, str(bed), out, motifs=[_motif("AAACCC")], p_value=1e-3
        )
        plus = [
            ln.split("\t")
            for ln in open(out).read().splitlines()
            if ln.endswith("+")
        ]
        assert plus and plus[0][1] == "103" and plus[0][2] == "109"

    def test_nested_output_dir_created(self, tmp_path):
        fasta = self._fasta(tmp_path, "TTTTTAAACCCTTTTT")
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t0\t16\n")
        out = str(tmp_path / "sub" / "dir" / "mpbs.bed")
        run_motif_matching(
            fasta, str(bed), out, motifs=[_motif("AAACCC")], p_value=1e-3
        )
        assert os.path.exists(out)
