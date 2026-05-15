"""Tests for deamtools.stat.plot_motif."""

import os

import numpy as np
import pandas as pd
import pysam
import pytest

from deamtools.stat.plot_motif import (
    _access_pwm,
    _atac_pwm,
    _pwm_to_information_df,
    _revcomp,
    run_plot_motif,
)

# Reference long enough to host motif windows comfortably (positions 0..29).
REF_SEQ = "ACGTACGTACGTACGTACGTACGTACGTAC"  # length 30
BAM_HEADER = {
    "HD": {"VN": "1.6"},
    "SQ": [{"LN": len(REF_SEQ), "SN": "chr1"}],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_read(name, seq, pos, *, is_reverse=False, mapq=30, baseq=40):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 0x10 if is_reverse else 0
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * len(seq))
    return a


def _write_bam(path, reads):
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=BAM_HEADER) as bam:
        for r in reads:
            bam.write(r)
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


# ---------------------------------------------------------------------------
# Pure-python helper tests
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_revcomp(self):
        assert _revcomp("ACGTN") == "NACGT"
        assert _revcomp("aacc") == "ggtt"

    def test_pwm_to_information_df_uniform_yields_zero_bits(self):
        # Equal counts at every position -> probability 0.25 each -> 0 bits.
        pwm = {b: [10.0] * 4 for b in "ACGTN"}
        df = _pwm_to_information_df(pwm, window_size=4)
        assert list(df.columns) == ["A", "C", "G", "T"]
        assert df.shape == (4, 4)
        assert np.allclose(df.values, 0.0)

    def test_pwm_to_information_df_pure_yields_two_bits(self):
        # All counts on A -> 2 bits at A, 0 elsewhere.
        pwm = {"A": [10.0] * 3, "C": [0.0] * 3, "G": [0.0] * 3,
               "T": [0.0] * 3, "N": [0.0] * 3}
        df = _pwm_to_information_df(pwm, window_size=3)
        assert np.allclose(df["A"], 2.0)
        assert np.allclose(df["C"], 0.0)
        assert np.allclose(df["G"], 0.0)
        assert np.allclose(df["T"], 0.0)
        # Index centred on zero
        assert list(df.index) == [-1, 0, 1]

    def test_pwm_to_information_df_zero_total_position(self):
        # A position with no observations -> all-zero row, no NaN.
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
# PWM builder tests
# ---------------------------------------------------------------------------


class TestAccessPWM:
    def test_forward_edit_excludes_centre(self, tmp_path, fasta_file):
        # REF_SEQ[5] == 'C'. Forward read of length 10 starting at pos 0,
        # with the read base at index 5 changed C -> T -> a deamination edit.
        assert REF_SEQ[5] == "C"
        seq = list(REF_SEQ[0:10])
        seq[5] = "T"
        read = _make_read("r1", "".join(seq), 0)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])

        regions = {"chr1": [(0, len(REF_SEQ))]}
        pwm = _access_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom=regions, window_size=4,
            min_mapq=0, min_baseq=0,
        )
        # Window size 4 centred at ref pos 5 -> p1=3, p2=7 -> REF_SEQ[3:7]="TACG".
        # Centre index = window_size // 2 = 2 (skipped).
        # Expected counts: pos 0 'T', pos 1 'A', pos 2 skipped, pos 3 'G'.
        assert pwm["T"][0] == 1 and sum(pwm[b][0] for b in "ACG") == 0
        assert pwm["A"][1] == 1 and sum(pwm[b][1] for b in "CGT") == 0
        assert all(pwm[b][2] == 0 for b in "ACGT")
        assert pwm["G"][3] == 1 and sum(pwm[b][3] for b in "ACT") == 0

    def test_no_edit_yields_empty_pwm(self, tmp_path, fasta_file):
        read = _make_read("r1", REF_SEQ[0:10], 0)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        pwm = _access_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
            window_size=4, min_mapq=0, min_baseq=0,
        )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")

    def test_reverse_edit_revcomps_window(self, tmp_path, fasta_file):
        # REF_SEQ[6] == 'G'. A reverse-strand read with read base 'A' at ref
        # pos 6 should register a G->A deamination event. The motif window
        # is then taken from the + strand and reverse-complemented.
        assert REF_SEQ[6] == "G"
        seq = list(REF_SEQ[0:10])
        seq[6] = "A"
        read = _make_read("r1", "".join(seq), 0, is_reverse=True)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])

        pwm = _access_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
            window_size=4, min_mapq=0, min_baseq=0,
        )
        # Window around pos 6, size 4 -> ref positions 4..7 = REF_SEQ[4:8] = "ACGT".
        # Reverse-complement: "ACGT" (palindromic). Centre = window_size//2 = 2.
        # Counts: pos 0 'A', pos 1 'C', pos 2 skipped, pos 3 'T'.
        assert pwm["A"][0] == 1 and sum(pwm[b][0] for b in "CGT") == 0
        assert pwm["C"][1] == 1 and sum(pwm[b][1] for b in "AGT") == 0
        assert all(pwm[b][2] == 0 for b in "ACGT")
        assert pwm["T"][3] == 1 and sum(pwm[b][3] for b in "ACG") == 0

    def test_baseq_filter_drops_edit(self, tmp_path, fasta_file):
        seq = list(REF_SEQ[0:10])
        seq[5] = "T"  # C->T at ref pos 5 (C)
        read = _make_read("r1", "".join(seq), 0, baseq=5)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        pwm = _access_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
            window_size=4, min_mapq=0, min_baseq=20,
        )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")

    def test_edit_too_close_to_chrom_end_skipped(self, tmp_path, fasta_file):
        # C->T edit at the very leftmost C (ref pos 1). Window of size 4 wants
        # positions -1..2; p1 < 0 so the edit must be skipped.
        assert REF_SEQ[1] == "C"
        seq = list(REF_SEQ[0:10])
        seq[1] = "T"
        read = _make_read("r1", "".join(seq), 0)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        pwm = _access_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
            window_size=4, min_mapq=0, min_baseq=0,
        )
        assert all(sum(pwm[b]) == 0 for b in "ACGT")


class TestATACPwm:
    def test_forward_cut_site_offset_by_4(self, tmp_path, fasta_file):
        # Forward read at pos 5, length 10. Cut site = 5 + 4 = 9.
        # window_size=4 -> p1 = 9 - 2 = 7, p2 = 11. REF_SEQ[7:11]="TACG".
        read = _make_read("r1", REF_SEQ[5:15], 5)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        pwm = _atac_pwm(
            bam_path=bam, fasta_path=fasta_file,
            regions_by_chrom={"chr1": [(0, len(REF_SEQ))]},
            window_size=4, min_mapq=0,
        )
        # Counts at positions 0..3: T, A, C, G  (no centre skip in atac mode).
        assert pwm["T"][0] == 1
        assert pwm["A"][1] == 1
        assert pwm["C"][2] == 1
        assert pwm["G"][3] == 1


# ---------------------------------------------------------------------------
# End-to-end run_plot_motif
# ---------------------------------------------------------------------------


class TestRunPlotMotif:
    def _seed_edit(self, tmp_path):
        # Build a BAM with several forward C->T edits at pos 4 to give the
        # whole-genome run something to plot.
        reads = []
        for i in range(5):
            seq = list(REF_SEQ[0:10])
            seq[4] = "T"
            reads.append(_make_read(f"r{i}", "".join(seq), 0))
        return _write_bam(str(tmp_path / "x.bam"), reads)

    def test_writes_plot_and_csv(self, tmp_path, fasta_file):
        bam = self._seed_edit(tmp_path)
        out = str(tmp_path / "motif.png")
        df = run_plot_motif(
            bam_path=bam, fasta_path=fasta_file, output_path=out,
            mode="access", window_size=4,
            min_mapq=0, min_baseq=0,
        )
        assert os.path.exists(out)
        assert os.path.exists(str(tmp_path / "motif.csv"))
        # Reload the CSV and verify shape
        loaded = pd.read_csv(tmp_path / "motif.csv", index_col=0)
        assert list(loaded.columns) == ["A", "C", "G", "T"]
        assert loaded.shape == df.shape == (4, 4)

    def test_invalid_mode_raises(self, tmp_path, fasta_file):
        bam = self._seed_edit(tmp_path)
        with pytest.raises(ValueError, match="mode"):
            run_plot_motif(
                bam_path=bam, fasta_path=fasta_file,
                output_path=str(tmp_path / "x.png"),
                mode="bogus", window_size=4,
            )

    def test_no_observations_raises(self, tmp_path, fasta_file):
        # A read with no edits -> access mode finds no centres -> error.
        read = _make_read("r1", REF_SEQ[0:10], 0)
        bam = _write_bam(str(tmp_path / "x.bam"), [read])
        with pytest.raises(RuntimeError, match="No editing"):
            run_plot_motif(
                bam_path=bam, fasta_path=fasta_file,
                output_path=str(tmp_path / "x.png"),
                mode="access", window_size=4,
                min_mapq=0, min_baseq=0,
            )

    def test_bed_restricts_processing(self, tmp_path, fasta_file):
        bam = self._seed_edit(tmp_path)
        bed = tmp_path / "r.bed"
        # Restrict to a region that does NOT contain pos 4 -> no edits found.
        bed.write_text("chr1\t10\t20\n")
        with pytest.raises(RuntimeError, match="No editing"):
            run_plot_motif(
                bam_path=bam, fasta_path=fasta_file,
                output_path=str(tmp_path / "x.png"),
                bed_path=str(bed), mode="access", window_size=4,
                min_mapq=0, min_baseq=0,
            )
