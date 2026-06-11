"""Tests for deamtools.footprint.footprint."""

import os

import pyBigWig
import pytest

from deamtools.footprint.footprint import _footprint_score, run_footprint


def _write_bigwig(path: str, chrom_len: int, entries: list[tuple[int, float]]) -> None:
    """Write a per-base BigWig; `entries` is a sorted list of (pos, value)."""
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", chrom_len)])
    positions = [p for p, _ in entries]
    values = [float(v) for _, v in entries]
    bw.addEntries("chr1", positions, values=values, span=1)
    bw.close()


def _rows(path: str) -> list[list[str]]:
    return [ln.split("\t") for ln in open(path).read().splitlines() if ln]


class TestFootprintScore:
    def test_score_formula(self):
        import numpy as np

        # L=2: left=[1,1], centre=[0,0], right=[1,1] -> 1 + 1 - 0 = 2
        signal = np.array([1.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        assert _footprint_score(signal, 2) == pytest.approx(2.0)


class TestRunFootprint:
    def test_footprint_detected(self, tmp_path):
        # Motif at [40,50) (L=10); window [30,60). Flanks=10, centre unwritten=0.
        bw = str(tmp_path / "s.bw")
        flank = [(p, 10.0) for p in range(30, 40)] + [(p, 10.0) for p in range(50, 60)]
        _write_bigwig(bw, 100, flank)
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\tTF1\n")

        run_footprint(bw, str(bed), str(tmp_path), "fp", n_shuffles=200, seed=0)

        rows = _rows(str(tmp_path / "fp.bed"))
        assert len(rows) == 1
        chrom, start, end, name, score, pval = rows[0]
        assert (chrom, start, end, name) == ("chr1", "40", "50", "TF1")
        assert float(score) == pytest.approx(20.0)   # 10 + 10 - 0
        assert float(pval) < 0.05                     # clear footprint

    def test_uniform_region_pvalue_one(self, tmp_path):
        # Uniform signal: score = level + level - level = level (here 5), but
        # every permutation gives the same score, so the p-value is 1.0.
        bw = str(tmp_path / "s.bw")
        _write_bigwig(bw, 100, [(p, 5.0) for p in range(0, 90)])
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\tTF1\n")

        run_footprint(bw, str(bed), str(tmp_path), "fp", n_shuffles=100, seed=0)
        chrom, start, end, name, score, pval = _rows(str(tmp_path / "fp.bed"))[0]
        assert float(score) == pytest.approx(5.0)
        assert float(pval) == pytest.approx(1.0)

    def test_negative_score_pvalue_one(self, tmp_path):
        # Centre high, flanks zero -> fp_score = 0 + 0 - 10 < 0 -> p_value 1.0
        # (the short-circuit, no permutations needed).
        bw = str(tmp_path / "s.bw")
        _write_bigwig(bw, 100, [(p, 10.0) for p in range(40, 50)])
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\tTF1\n")

        run_footprint(bw, str(bed), str(tmp_path), "fp", n_shuffles=100, seed=0)
        chrom, start, end, name, score, pval = _rows(str(tmp_path / "fp.bed"))[0]
        assert float(score) == pytest.approx(-10.0)
        assert float(pval) == pytest.approx(1.0)

    def test_out_of_bounds_region_skipped(self, tmp_path):
        bw = str(tmp_path / "s.bw")
        _write_bigwig(bw, 100, [(p, 10.0) for p in range(0, 100)])
        bed = tmp_path / "r.bed"
        # window for [95,99) (L=4) is [91,103) -> exceeds chrom length 100.
        bed.write_text("chr1\t95\t99\tTF1\n")

        run_footprint(bw, str(bed), str(tmp_path), "fp", n_shuffles=10, seed=0)
        assert _rows(str(tmp_path / "fp.bed")) == []

    def test_name_defaults_to_dot_without_column4(self, tmp_path):
        bw = str(tmp_path / "s.bw")
        flank = [(p, 8.0) for p in range(30, 40)] + [(p, 8.0) for p in range(50, 60)]
        _write_bigwig(bw, 100, flank)
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\n")  # no name column

        run_footprint(bw, str(bed), str(tmp_path), "fp", n_shuffles=50, seed=0)
        assert _rows(str(tmp_path / "fp.bed"))[0][3] == "."

    def test_seed_makes_pvalues_reproducible(self, tmp_path):
        bw = str(tmp_path / "s.bw")
        flank = [(p, 3.0) for p in range(30, 40)] + [(p, 3.0) for p in range(50, 60)]
        _write_bigwig(bw, 100, flank)
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\tTF1\n")

        run_footprint(bw, str(bed), str(tmp_path), "a", n_shuffles=200, seed=42)
        run_footprint(bw, str(bed), str(tmp_path), "b", n_shuffles=200, seed=42)
        assert _rows(str(tmp_path / "a.bed")) == _rows(str(tmp_path / "b.bed"))

    def test_nested_out_dir_created(self, tmp_path):
        bw = str(tmp_path / "s.bw")
        _write_bigwig(bw, 100, [(p, 1.0) for p in range(0, 90)])
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t40\t50\tTF1\n")
        out_dir = str(tmp_path / "sub" / "dir")
        run_footprint(bw, str(bed), out_dir, "fp", n_shuffles=10, seed=0)
        assert os.path.exists(os.path.join(out_dir, "fp.bed"))
