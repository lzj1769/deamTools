"""Tests for deamtools.stat.plot_motif (BigWig-driven)."""

import os

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import pytest

from deamtools.stat.plot_motif import (
    _accumulate_pwm_from_bigwig,
    _pwm_to_information_df,
    _revcomp,
    run_plot_motif,
)

# Reference long enough to host motif windows comfortably (positions 0..29).
# Pattern: A C G T A C G T ...  -- C at every odd index, G at every (i % 4 == 2).
REF_SEQ = "ACGTACGTACGTACGTACGTACGTACGTAC"  # length 30


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_fasta(tmp_path, seq=REF_SEQ, chrom="chr1"):
    path = str(tmp_path / "ref.fa")
    with open(path, "w") as f:
        f.write(f">{chrom}\n{seq}\n")
    pysam.faidx(path)
    return path


def _write_bigwig(tmp_path, entries, chrom_sizes=None):
    """Write a BigWig file at single-base resolution.

    Parameters
    ----------
    entries : list[tuple[str, int, float]]
        Per-base entries to write, given as ``(chrom, position, value)``
        with value > 0. The chromosome must appear in ``chrom_sizes``.
    chrom_sizes : dict[str, int], optional
        Defaults to ``{"chr1": len(REF_SEQ)}``.
    """
    if chrom_sizes is None:
        chrom_sizes = {"chr1": len(REF_SEQ)}
    path = str(tmp_path / "s.bw")
    bw = pyBigWig.open(path, "w")
    bw.addHeader(list(chrom_sizes.items()))
    # pyBigWig requires entries to be added in (chrom, start) order.
    entries = sorted(entries, key=lambda t: (list(chrom_sizes).index(t[0]), t[1]))
    for chrom, pos, value in entries:
        bw.addEntries(chrom, [pos], values=[float(value)], span=1)
    bw.close()
    return path


# ---------------------------------------------------------------------------
# Pure-python helper tests
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_revcomp(self):
        assert _revcomp("ACGTN") == "NACGT"
        assert _revcomp("aacc") == "ggtt"

    def test_pwm_to_information_df_uniform_yields_zero_bits(self):
        pwm = {b: [10.0] * 4 for b in "ACGTN"}
        df = _pwm_to_information_df(pwm, window_size=4)
        assert list(df.columns) == ["A", "C", "G", "T"]
        assert df.shape == (4, 4)
        assert np.allclose(df.values, 0.0)

    def test_pwm_to_information_df_pure_yields_two_bits(self):
        pwm = {"A": [10.0] * 3, "C": [0.0] * 3, "G": [0.0] * 3,
               "T": [0.0] * 3, "N": [0.0] * 3}
        df = _pwm_to_information_df(pwm, window_size=3)
        assert np.allclose(df["A"], 2.0)
        assert np.allclose(df["C"], 0.0)
        assert np.allclose(df["G"], 0.0)
        assert np.allclose(df["T"], 0.0)
        assert list(df.index) == [-1, 0, 1]

    def test_pwm_to_information_df_zero_total_position(self):
        pwm = {b: [0.0, 5.0] for b in "ACGTN"}
        pwm["A"][1] = 5.0
        pwm["C"][1] = 0.0
        pwm["G"][1] = 0.0
        pwm["T"][1] = 0.0
        pwm["N"][1] = 0.0
        df = _pwm_to_information_df(pwm, window_size=2)
        assert not df.isna().any().any()
        assert np.allclose(df.iloc[0].values, 0.0)


# ---------------------------------------------------------------------------
# _accumulate_pwm_from_bigwig
# ---------------------------------------------------------------------------


class TestAccumulatePWM:
    def test_forward_edit_at_C_excludes_centre(self, tmp_path):
        # REF_SEQ[5] == 'C'. A count of 1 at position 5 should contribute the
        # surrounding window (positions 3..6 for window_size=4) to the PWM with
        # the centre column (index 2) left at zero.
        assert REF_SEQ[5] == "C"
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 5, 1.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
                window_size=4,
            )

        # Window size 4 centred at pos 5 -> REF_SEQ[3:7] = "TACG".
        # Centre index = window_size // 2 = 2 (skipped).
        # Expected: pos 0 'T', pos 1 'A', pos 2 skipped, pos 3 'G'.
        assert pwm["T"][0] == 1 and sum(pwm[b][0] for b in "ACG") == 0
        assert pwm["A"][1] == 1 and sum(pwm[b][1] for b in "CGT") == 0
        assert all(pwm[b][2] == 0 for b in "ACGT")
        assert pwm["G"][3] == 1 and sum(pwm[b][3] for b in "ACT") == 0

    def test_reverse_edit_at_G_revcomps_window(self, tmp_path):
        # REF_SEQ[6] == 'G'. A count at position 6 should be treated as a
        # reverse-strand deamination event: window taken from + strand and
        # reverse-complemented before being accumulated.
        assert REF_SEQ[6] == "G"
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 6, 1.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
                window_size=4,
            )

        # Window around pos 6, size 4 -> REF_SEQ[4:8] = "ACGT".
        # Reverse-complement: "ACGT" (palindromic in this stretch).
        # Centre index = 2 (skipped).
        # Counts: pos 0 'A', pos 1 'C', pos 2 skipped, pos 3 'T'.
        assert pwm["A"][0] == 1 and sum(pwm[b][0] for b in "CGT") == 0
        assert pwm["C"][1] == 1 and sum(pwm[b][1] for b in "AGT") == 0
        assert all(pwm[b][2] == 0 for b in "ACGT")
        assert pwm["T"][3] == 1 and sum(pwm[b][3] for b in "ACG") == 0

    def test_count_acts_as_weight(self, tmp_path):
        # A count of 7 at a single C should contribute 7 to each flanking
        # column (mathematically identical to seven independent unit events
        # at the same position).
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 5, 7.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
                window_size=4,
            )
        # Same flanking pattern as test_forward_edit_at_C_excludes_centre,
        # but with weight 7.
        assert pwm["T"][0] == 7
        assert pwm["A"][1] == 7
        assert all(pwm[b][2] == 0 for b in "ACGT")
        assert pwm["G"][3] == 7

    def test_non_CG_centre_skipped(self, tmp_path):
        # A count at a reference 'A' position must not contribute (it cannot
        # be a deamination event). REF_SEQ[0] == 'A'.
        assert REF_SEQ[0] == "A"
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 0, 5.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
                window_size=4,
            )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")

    def test_edit_too_close_to_chrom_end_skipped(self, tmp_path):
        # REF_SEQ[1] == 'C'. Window of size 4 centred at pos 1 wants positions
        # -1..2; p1 < 0 so this site must be skipped.
        assert REF_SEQ[1] == "C"
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 1, 1.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
                window_size=4,
            )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")

    def test_region_restriction(self, tmp_path):
        # Two count-1 entries, only one within the requested region. The
        # in-region site should contribute; the out-of-region one should not.
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path,
                                 entries=[("chr1", 5, 1.0), ("chr1", 13, 1.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr1": [(0, 10)]},  # excludes pos 13
                window_size=4,
            )
        # Only one site contributed -> totals at each column should be 1.
        assert pwm["T"][0] == 1
        assert pwm["A"][1] == 1
        assert pwm["G"][3] == 1

    def test_missing_chromosome_skipped(self, tmp_path):
        # A region on a chromosome that is not present in the BigWig header
        # should be silently skipped (with a warning).
        fasta = _write_fasta(tmp_path)
        bw_path = _write_bigwig(tmp_path, entries=[("chr1", 5, 1.0)])

        with pyBigWig.open(bw_path) as bw, pysam.FastaFile(fasta) as fa:
            pwm = _accumulate_pwm_from_bigwig(
                bw=bw, fasta=fa,
                regions_by_chrom={"chr_missing": [(0, 100)]},
                window_size=4,
            )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")


# ---------------------------------------------------------------------------
# End-to-end run_plot_motif
# ---------------------------------------------------------------------------


class TestRunPlotMotif:
    def _seed_bigwig(self, tmp_path):
        # Five count-1 entries at REF_SEQ[5] == 'C'. With count weighting,
        # this is equivalent to a single entry of weight 5 (we use weight 5
        # directly).
        return _write_bigwig(tmp_path, entries=[("chr1", 5, 5.0)])

    def test_writes_plot_and_csv(self, tmp_path):
        fasta = _write_fasta(tmp_path)
        bw = self._seed_bigwig(tmp_path)
        out = str(tmp_path / "motif.png")
        df = run_plot_motif(
            bigwig_path=bw, fasta_path=fasta, output_path=out,
            window_size=4,
        )
        assert os.path.exists(out)
        assert os.path.exists(str(tmp_path / "motif.csv"))
        loaded = pd.read_csv(tmp_path / "motif.csv", index_col=0)
        assert list(loaded.columns) == ["A", "C", "G", "T"]
        assert loaded.shape == df.shape == (4, 4)

    def test_invalid_window_size_raises(self, tmp_path):
        fasta = _write_fasta(tmp_path)
        bw = self._seed_bigwig(tmp_path)
        with pytest.raises(ValueError, match="window_size"):
            run_plot_motif(
                bigwig_path=bw, fasta_path=fasta,
                output_path=str(tmp_path / "x.png"),
                window_size=1,
            )

    def test_no_observations_raises(self, tmp_path):
        # Empty BigWig (no entries) -> no editing sites -> error.
        fasta = _write_fasta(tmp_path)
        bw = _write_bigwig(tmp_path, entries=[])
        with pytest.raises(RuntimeError, match="No editing"):
            run_plot_motif(
                bigwig_path=bw, fasta_path=fasta,
                output_path=str(tmp_path / "x.png"),
                window_size=4,
            )

    def test_bed_restricts_processing(self, tmp_path):
        fasta = _write_fasta(tmp_path)
        bw = self._seed_bigwig(tmp_path)
        bed = tmp_path / "r.bed"
        # Restrict to a region that does NOT contain pos 5 -> no edits found.
        bed.write_text("chr1\t10\t20\n")
        with pytest.raises(RuntimeError, match="No editing"):
            run_plot_motif(
                bigwig_path=bw, fasta_path=fasta,
                output_path=str(tmp_path / "x.png"),
                bed_path=str(bed),
                window_size=4,
            )
