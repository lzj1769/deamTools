"""Tests for deamtools.preprocessing.bam2bw."""

import os

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import pytest

from deamtools.preprocessing.bam2bw import (
    _signal_for_region,
    run_bam2bw,
)
from deamtools.utils.regions import _load_regions

# ---------------------------------------------------------------------------
# Reference sequence used across all tests.
# C positions (forward): 1, 4, 8
# G positions (reverse):  2, 5, 9
# ---------------------------------------------------------------------------
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
    ref_id: int,
    pos: int,
    mapq: int = 30,
    is_reverse: bool = False,
    baseq: int = 40,
    extra_flags: int = 0,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = extra_flags | (0x10 if is_reverse else 0)
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * len(seq))
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


@pytest.fixture()
def chrom_sizes_file(tmp_path):
    path = tmp_path / "chrom.sizes"
    path.write_text(f"chr1\t{len(REF_SEQ)}\n")
    return str(path)


# ---------------------------------------------------------------------------
# _load_regions
# ---------------------------------------------------------------------------


def _expect(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


class TestLoadRegions:
    def test_basic(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\n")
        result = _load_regions(str(bed))
        assert isinstance(result, pd.DataFrame)
        pd.testing.assert_frame_equal(result, _expect([("chr1", 0, 5)]))

    def test_overlapping_intervals_merged(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\nchr1\t3\t8\n")
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)), _expect([("chr1", 0, 8)])
        )

    def test_adjacent_intervals_merged(self, tmp_path):
        # End of one == start of next: BED is half-open, but treating these as
        # a contiguous region prevents double-counting reads at the boundary.
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\nchr1\t5\t8\n")
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)), _expect([("chr1", 0, 8)])
        )

    def test_non_overlapping_intervals_kept(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t3\nchr1\t5\t8\n")
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)),
            _expect([("chr1", 0, 3), ("chr1", 5, 8)]),
        )

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("# header\nchr1\t0\t5\n")
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)), _expect([("chr1", 0, 5)])
        )

    def test_track_and_browser_lines_skipped(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text(
            "browser position chr1:1-1000\n"
            'track name="x" description="y"\n'
            "chr1\t0\t5\n"
        )
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)), _expect([("chr1", 0, 5)])
        )

    def test_extra_columns_ignored(self, tmp_path):
        # A 6-column BED line (chrom, start, end, name, score, strand) loads.
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t10\t20\tregion_a\t900\t+\n")
        pd.testing.assert_frame_equal(
            _load_regions(str(bed)), _expect([("chr1", 10, 20)])
        )

    def test_multiple_chromosomes_sorted(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr2\t10\t20\nchr1\t0\t5\n")
        result = _load_regions(str(bed))
        assert set(result["chrom"]) == {"chr1", "chr2"}
        # Rows should be sorted by chromosome then start.
        assert list(result["chrom"]) == ["chr1", "chr2"]

    def test_empty_file_returns_empty_dataframe(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("# only headers\n\n")
        result = _load_regions(str(bed))
        assert list(result.columns) == ["chrom", "start", "end"]
        assert len(result) == 0

    def test_invalid_interval_start_gt_end_raises(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t10\t5\n")
        with pytest.raises(ValueError, match="start > end"):
            _load_regions(str(bed))


# ---------------------------------------------------------------------------
# _count_deamination_on_chrom
# ---------------------------------------------------------------------------


def _run_count_region(
    bam_path, fasta_file, regions, *, min_mapq, min_baseq, extend_size
):
    """Stitch per-region count signal into a full-chromosome array.

    Lets the legacy tests in :class:`TestCountDeamination` keep asserting on
    absolute reference positions even though :func:`_signal_for_region`
    returns a per-region array.
    """
    if regions is None:
        regions = [(0, len(REF_SEQ))]
    out = [0.0] * len(REF_SEQ)
    for start, end in regions:
        _, _, _, sig = _signal_for_region(
            bam_path=bam_path,
            fasta_path=fasta_file,
            chrom="chr1",
            start=start,
            end=end,
            mode="count",
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
            min_coverage=0,
        )
        for i, v in enumerate(sig):
            out[start + i] = float(v)
    import numpy as _np

    return _np.array(out, dtype=_np.float32)


class TestCountDeamination:
    """Unit tests for the per-region count signal."""

    def _run(
        self,
        bam_path,
        fasta_file,
        *,
        regions=None,
        min_mapq=0,
        min_baseq=0,
        extend_size=0,
    ):
        return _run_count_region(
            bam_path,
            fasta_file,
            regions,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
        )

    # -- Deamination detection ------------------------------------------------

    def test_forward_ct_event_detected(self, tmp_path, fasta_file):
        # REF:  A C G T C G A T C G
        # READ: A T G T C G A T C G  ← C→T at pos 1
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file)
        assert counts[1] == 1.0
        assert counts.sum() == 1.0

    def test_reverse_ga_event_detected(self, tmp_path, fasta_file):
        # Reverse-strand read; G→A at ref pos 5 (- strand C deaminated to T).
        # query_sequence derived by reverse-complementing the deaminated - strand:
        #   normal rev-comp of REF = "CGATCGACGT"
        #   with deamination at ref pos 5: raw read = "CGATTGACGT"
        #   stored in BAM (rev-comp of raw) = "ACGTCAATCG"
        bam = _write_bam(
            str(tmp_path / "t.bam"),
            [_make_read("r1", "ACGTCAATCG", 0, 0, is_reverse=True)],
        )
        counts = self._run(bam, fasta_file)
        assert counts[5] == 1.0
        assert counts.sum() == 1.0

    def test_no_event_on_exact_match(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"), [_make_read("r1", REF_SEQ, 0, 0)])
        assert self._run(bam, fasta_file).sum() == 0.0

    def test_multiple_reads_accumulate(self, tmp_path, fasta_file):
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),
            _make_read("r2", "ATGTCGATCG", 0, 0),
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        assert self._run(bam, fasta_file)[1] == 2.0

    def test_multiple_events_in_one_read(self, tmp_path, fasta_file):
        # REF:  A C G T C G A T C G
        # READ: A T G T T G A T C G  ← C→T at pos 1 and pos 4
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTTGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file)
        assert counts[1] == 1.0
        assert counts[4] == 1.0
        assert counts.sum() == 2.0

    # -- Filtering ------------------------------------------------------------

    def test_low_mapq_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0, mapq=10)]
        )
        assert self._run(bam, fasta_file, min_mapq=20).sum() == 0.0

    def test_mapq_at_threshold_included(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0, mapq=20)]
        )
        assert self._run(bam, fasta_file, min_mapq=20).sum() == 1.0

    def test_low_baseq_position_excluded(self, tmp_path, fasta_file):
        # baseq=5 < min_baseq=20 → the C→T event at pos 1 is not counted
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0, baseq=5)]
        )
        assert self._run(bam, fasta_file, min_baseq=20).sum() == 0.0

    def test_secondary_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"),
            [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x100)],
        )
        assert self._run(bam, fasta_file).sum() == 0.0

    def test_duplicate_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"),
            [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x400)],
        )
        assert self._run(bam, fasta_file).sum() == 0.0

    def test_supplementary_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"),
            [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x800)],
        )
        assert self._run(bam, fasta_file).sum() == 0.0

    # -- Region restriction ---------------------------------------------------

    def test_event_inside_region_counted(self, tmp_path, fasta_file):
        # C→T at pos 1 and pos 4; restrict to [0, 3) → only pos 1 counted
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTTGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file, regions=[(0, 3)])
        assert counts[1] == 1.0
        assert counts[4] == 0.0
        assert counts.sum() == 1.0

    def test_event_outside_region_excluded(self, tmp_path, fasta_file):
        # C→T at pos 1; restrict to [5, 10) → no events in region
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file, regions=[(5, 10)])
        assert counts.sum() == 0.0

    # -- Signal extension -----------------------------------------------------

    def test_extend_size_spreads_signal(self, tmp_path, fasta_file):
        # Single C→T at pos 1; extend_size=2 → signal at positions 0,1,2,3
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file, extend_size=2)
        assert counts[0] == 1.0
        assert counts[1] == 1.0
        assert counts[2] == 1.0
        assert counts[3] == 1.0
        assert counts[4] == 0.0

    def test_extend_size_zero_unchanged(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        counts = self._run(bam, fasta_file, extend_size=0)
        assert counts[1] == 1.0
        assert counts.sum() == 1.0


# ---------------------------------------------------------------------------
# run_bam2bw (integration)
# ---------------------------------------------------------------------------


class TestRunBam2bw:
    """End-to-end tests that write a real BigWig and check its content."""

    def test_bigwig_created_with_correct_signal(
        self, tmp_path, fasta_file, chrom_sizes_file
    ):
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
        )

        assert os.path.exists(out)
        with pyBigWig.open(out) as bw:
            # Event at pos 1
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(1.0)
            # No event at pos 0
            assert bw.stats("chr1", 0, 1, type="mean")[0] is None

    def test_normalize_count_scales_by_total(
        self, tmp_path, fasta_file, chrom_sizes_file
    ):
        # Read with C->T at pos 1 and pos 4 -> 2 edits genome-wide.
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTTGATCG", 0, 0)]
        )
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
            normalize=True,
            scale_factor=10.0,
        )

        with pyBigWig.open(out) as bw:
            # raw count 1 -> 1 * scale_factor / total = 1 * 10 / 2 = 5.0
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(5.0)
            assert bw.stats("chr1", 4, 5, type="mean")[0] == pytest.approx(5.0)

    def test_normalize_ignored_in_ratio_mode(
        self, tmp_path, fasta_file, chrom_sizes_file
    ):
        # 1 edit (T) + 1 ref C at pos 1 -> ratio 0.5; --normalize must not change it.
        reads = [_make_read("r1", "ATGTCGATCG", 0, 0), _make_read("r2", REF_SEQ, 0, 0)]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
            mode="ratio",
            min_coverage=0,
            normalize=True,
            scale_factor=10.0,
        )
        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(0.5)

    def test_infer_chrom_sizes_from_bam(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=None,
            min_mapq=0,
            min_baseq=0,
        )
        assert os.path.exists(out)

    def test_bed_restricts_output_signal(self, tmp_path, fasta_file, chrom_sizes_file):
        # C→T at pos 1 and pos 4; BED restricts to [0, 3)
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTTGATCG", 0, 0)]
        )
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t3\n")
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            bed_path=str(bed),
            min_mapq=0,
            min_baseq=0,
        )

        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(1.0)
            assert bw.stats("chr1", 4, 5, type="mean")[0] is None

    def test_output_parent_directory_created(
        self, tmp_path, fasta_file, chrom_sizes_file
    ):
        bam = _write_bam(str(tmp_path / "t.bam"), [_make_read("r1", REF_SEQ, 0, 0)])
        out = str(tmp_path / "nested" / "dir" / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path / "nested" / "dir"),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
        )
        assert os.path.exists(out)

    def test_no_events_writes_empty_bigwig(
        self, tmp_path, fasta_file, chrom_sizes_file
    ):
        # Read identical to reference → no deamination events
        bam = _write_bam(str(tmp_path / "t.bam"), [_make_read("r1", REF_SEQ, 0, 0)])
        out = str(tmp_path / "out.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="out",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
        )

        assert os.path.exists(out)
        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 0, len(REF_SEQ), type="mean")[0] is None


# ---------------------------------------------------------------------------
# Ratio mode
# ---------------------------------------------------------------------------


class TestRatioMode:
    """Tests for --mode ratio: edit_count / total ACGT coverage.

    Mirrors the ACCESS-ATAC reference algorithm: the denominator is the
    sum of ACGT read counts at each position (not only "informative"
    bases), and positions below ``--min_coverage`` are masked to 0.
    """

    def _run(self, bam_path, fasta_file, *, extend_size=0, min_coverage=0):
        _, _, _, signal = _signal_for_region(
            bam_path=bam_path,
            fasta_path=fasta_file,
            chrom="chr1",
            start=0,
            end=len(REF_SEQ),
            mode="ratio",
            min_mapq=0,
            min_baseq=0,
            extend_size=extend_size,
            min_coverage=min_coverage,
        )
        return signal

    def test_invalid_mode_raises(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"), [_make_read("r1", REF_SEQ, 0, 0)])
        with pytest.raises(ValueError, match="mode"):
            run_bam2bw(
                bam_path=bam,
                fasta_path=fasta_file,
                out_dir=str(tmp_path),
                out_name="out",
                min_mapq=0,
                min_baseq=0,
                mode="bogus",
            )

    def test_ratio_three_quarters_at_C(self, tmp_path, fasta_file):
        # Reference C at pos 1: 3 reads with T (events), 1 read with C.
        # Total ACGT coverage at pos 1 = 4. Ratio = 3 / 4 = 0.75.
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),
            _make_read("r2", "ATGTCGATCG", 0, 0),
            _make_read("r3", "ATGTCGATCG", 0, 0),
            _make_read("r4", REF_SEQ, 0, 0),
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file)
        assert ratio[1] == pytest.approx(0.75)
        # Position 4 (ref C, no edits, but 4 reads of coverage): 0/4 = 0.
        assert ratio[4] == 0.0

    def test_ratio_uses_total_coverage_denominator(self, tmp_path, fasta_file):
        # At ref C (pos 1): 1 read T (event), 1 read A (sequencing error/SNP).
        # Total ACGT coverage = 2 -> ratio = 1/2 = 0.5 under the reference
        # algorithm (not 1.0 as under an informative-only denominator).
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at C -> event
            _make_read("r2", "AAGTCGATCG", 0, 0),  # A at C -> coverage only
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file)
        assert ratio[1] == pytest.approx(0.5)

    def test_ratio_zero_when_no_edits(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"), [_make_read("r1", REF_SEQ, 0, 0)])
        ratio = self._run(bam, fasta_file)
        assert (ratio == 0.0).all()

    def test_ratio_reverse_strand_strand_agnostic(self, tmp_path, fasta_file):
        # Reverse-strand reads at ref G (pos 2). Two reads with A at ref pos 2
        # contribute events; one with G contributes coverage only.
        # Total coverage at pos 2 = 3, events = 2, ratio = 2/3.
        seq = "ACATCGATCG"  # A at pos 2 (G->A reference mismatch)
        reads = [
            _make_read("r1", seq, 0, 0, is_reverse=True),
            _make_read("r2", seq, 0, 0, is_reverse=True),
            _make_read("r3", REF_SEQ, 0, 0, is_reverse=True),
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file)
        assert ratio[2] == pytest.approx(2.0 / 3.0)

    def test_ratio_extend_size_ignored(self, tmp_path, fasta_file):
        # Per the reference algorithm, --extend_size only applies to count
        # mode. With a single C->T at pos 1 and extend_size=2, the ratio
        # mode signal must still be non-zero only at pos 1 (the editing site).
        bam = _write_bam(
            str(tmp_path / "t.bam"), [_make_read("r1", "ATGTCGATCG", 0, 0)]
        )
        ratio = self._run(bam, fasta_file, extend_size=2)
        assert ratio[1] == pytest.approx(1.0)
        # Neighbouring bases stay at zero -- no convolution / extension.
        assert ratio[0] == 0.0
        assert ratio[2] == 0.0
        assert ratio[3] == 0.0

    def test_min_coverage_threshold_masks_low_coverage_positions(
        self, tmp_path, fasta_file
    ):
        # 3 reads, all with C->T at pos 1. Total coverage = 3.
        reads = [_make_read(f"r{i}", "ATGTCGATCG", 0, 0) for i in range(3)]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        # Coverage threshold above actual coverage -> ratio masked to 0.
        ratio_masked = self._run(bam, fasta_file, min_coverage=10)
        assert ratio_masked[1] == 0.0
        # Coverage threshold at or below actual coverage -> ratio reported.
        ratio_kept = self._run(bam, fasta_file, min_coverage=3)
        assert ratio_kept[1] == pytest.approx(1.0)

    def test_ratio_bigwig_round_trip(self, tmp_path, fasta_file, chrom_sizes_file):
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at C (event)
            _make_read("r2", "ATGTCGATCG", 0, 0),  # T at C (event)
            _make_read("r3", REF_SEQ, 0, 0),  # C at C (coverage only)
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        out = str(tmp_path / "ratio.bw")
        run_bam2bw(
            bam_path=bam,
            fasta_path=fasta_file,
            out_dir=str(tmp_path),
            out_name="ratio",
            chrom_sizes_path=chrom_sizes_file,
            min_mapq=0,
            min_baseq=0,
            mode="ratio",
            min_coverage=0,
        )
        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(2.0 / 3.0)


# ---------------------------------------------------------------------------
# Fragment merging (paired-end)
# ---------------------------------------------------------------------------


def _make_pair(name, seq1, seq2, pos1=0, pos2=0, mapq=30, baseq=40):
    """A proper pair whose two mates can be merged into one fragment."""
    r1 = _make_read(
        name, seq1, 0, pos1, mapq=mapq, baseq=baseq, extra_flags=0x1 | 0x2 | 0x40
    )
    r1.next_reference_id = 0
    r1.next_reference_start = pos2
    r2 = _make_read(
        name,
        seq2,
        0,
        pos2,
        mapq=mapq,
        baseq=baseq,
        is_reverse=True,
        extra_flags=0x1 | 0x2 | 0x80,
    )
    r2.next_reference_id = 0
    r2.next_reference_start = pos1
    return [r1, r2]


class TestFragmentMerging:
    """A position both mates cover is one event, not two.

    REF_SEQ is ACGTCGATCG; position 1 is a reference C. Both mates of a real
    pair report the same edit there, because library prep turns the deaminated
    C into a T:A pair that both strands carry.
    """

    EDITED = "ATGTCGATCG"  # C->T at pos 1

    def _signal(self, bam_path, fasta_file, *, mode="count", min_mapq=0):
        _, _, _, signal = _signal_for_region(
            bam_path=bam_path,
            fasta_path=fasta_file,
            chrom="chr1",
            start=0,
            end=len(REF_SEQ),
            mode=mode,
            min_mapq=min_mapq,
            min_baseq=0,
            extend_size=0,
            min_coverage=0,
        )
        return signal

    def test_overlapping_mates_count_one_event(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "p.bam"), _make_pair("frag", self.EDITED, self.EDITED)
        )
        # Per record this was 2; the mates cover the same position, so it is 1.
        assert self._signal(bam, fasta_file)[1] == pytest.approx(1.0)

    def test_non_overlapping_mates_both_count(self, tmp_path, fasta_file):
        # R1 covers 0-4 with the C->T at 1, R2 covers 5-9 with the G->A at 5.
        reads = _make_pair("frag", "ATGTC", "AATCG", pos1=0, pos2=5)
        bam = _write_bam(str(tmp_path / "p.bam"), reads)
        signal = self._signal(bam, fasta_file)
        assert signal[1] == pytest.approx(1.0)
        assert signal[5] == pytest.approx(1.0)

    def test_extend_size_broadcasts_the_merged_event_once(self, tmp_path, fasta_file):
        bam = _write_bam(
            str(tmp_path / "p.bam"), _make_pair("frag", self.EDITED, self.EDITED)
        )
        _, _, _, signal = _signal_for_region(
            bam_path=bam,
            fasta_path=fasta_file,
            chrom="chr1",
            start=0,
            end=len(REF_SEQ),
            mode="count",
            min_mapq=0,
            min_baseq=0,
            extend_size=1,
            min_coverage=0,
        )
        # One event spread over pos 0-2, at height 1 rather than 2.
        assert signal[0] == pytest.approx(1.0)
        assert signal[1] == pytest.approx(1.0)
        assert signal[2] == pytest.approx(1.0)
        assert signal[3] == 0.0

    def test_ratio_denominator_is_also_per_fragment(self, tmp_path, fasta_file):
        # One edited fragment (two overlapping mates) plus one unedited single
        # read. Per record that was 2 edits over 3 reads = 0.667; per fragment
        # it is 1 edit over 2 fragments = 0.5.
        reads = _make_pair("frag", self.EDITED, self.EDITED)
        reads.append(_make_read("solo", REF_SEQ, 0, 0))
        bam = _write_bam(str(tmp_path / "p.bam"), reads)
        assert self._signal(bam, fasta_file, mode="ratio")[1] == pytest.approx(0.5)

    def test_ratio_denominator_honours_min_mapq(self, tmp_path, fasta_file):
        # The old denominator came from count_coverage, which applies no MAPQ
        # filter, so this low-MAPQ read used to dilute the ratio to 0.5 while
        # being excluded from the numerator.
        reads = [
            _make_read("keep", self.EDITED, 0, 0, mapq=30),
            _make_read("drop", REF_SEQ, 0, 0, mapq=0),
        ]
        bam = _write_bam(str(tmp_path / "p.bam"), reads)
        ratio = self._signal(bam, fasta_file, mode="ratio", min_mapq=20)
        assert ratio[1] == pytest.approx(1.0)

    def test_orphaned_mate_still_counts(self, tmp_path, fasta_file):
        # R2 fails --min_mapq, so R1's partner never arrives; R1 must still be
        # counted rather than sitting in the buffer.
        reads = _make_pair("frag", self.EDITED, self.EDITED)
        reads[1].mapping_quality = 0
        bam = _write_bam(str(tmp_path / "p.bam"), reads)
        assert self._signal(bam, fasta_file, min_mapq=20)[1] == pytest.approx(1.0)

    def test_single_end_reads_are_unaffected(self, tmp_path, fasta_file):
        reads = [_make_read(f"r{i}", self.EDITED, 0, 0) for i in range(3)]
        bam = _write_bam(str(tmp_path / "s.bam"), reads)
        assert self._signal(bam, fasta_file)[1] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# Tn5 cut sites (--event tn5)
# ---------------------------------------------------------------------------


class TestTn5Cuts:
    """--event tn5 counts insertion sites: each read's 5' end, shifted +4/-5.

    The shift puts both reads of one insertion on the same base. Tn5 nicks the
    two strands 9 bp apart, so the fragments either side of an insertion share a
    9-bp duplication [p, p + 9): a forward read of the right-hand fragment starts
    at p, a reverse read of the left-hand one ends (exclusive) at p + 9, and both
    must land on p + 4.
    """

    LENGTH = 200
    READ = 30

    def _bam(self, tmp_path, reads, name="t.bam"):
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": self.LENGTH}]}
        path = str(tmp_path / name)
        tmp = path + ".u.bam"
        with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
            for r in reads:
                bam.write(r)
        pysam.sort("-o", path, tmp)
        os.remove(tmp)
        pysam.index(path)
        return path

    def _fasta(self, tmp_path):
        path = str(tmp_path / "long.fa")
        with open(path, "w") as f:
            f.write(">chr1\n" + "ACGT" * (self.LENGTH // 4) + "\n")
        pysam.faidx(path)
        return path

    def _read(
        self,
        name,
        pos,
        *,
        reverse=False,
        mapq=30,
        baseq=40,
        flags=0,
        mate=None,
        read1=True,
    ):
        seq = "A" * self.READ
        a = _make_read(
            name,
            seq,
            0,
            pos,
            mapq=mapq,
            is_reverse=reverse,
            baseq=baseq,
            extra_flags=flags,
        )
        if mate is not None:
            a.flag |= 0x1 | 0x2 | (0x40 if read1 else 0x80)
            a.next_reference_id, a.next_reference_start = 0, mate
        return a

    def _cuts(
        self,
        tmp_path,
        reads,
        *,
        start=0,
        end=None,
        extend_size=0,
        min_mapq=0,
        min_baseq=0,
    ):
        bam = self._bam(tmp_path, reads)
        _, _, _, signal = _signal_for_region(
            bam_path=bam,
            fasta_path=self._fasta(tmp_path),
            chrom="chr1",
            start=start,
            end=self.LENGTH if end is None else end,
            mode="count",
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
            min_coverage=0,
            event="tn5",
        )
        return {start + int(i): float(signal[i]) for i in np.nonzero(signal)[0]}

    def test_both_sides_of_one_insertion_land_on_the_same_base(self, tmp_path):
        p = 100  # the 9-bp duplication is [100, 109)
        right = self._read("right", p)  # forward, starts at p
        left = self._read("left", p + 9 - self.READ, reverse=True)  # ends at p+9
        assert left.reference_end == p + 9
        assert self._cuts(tmp_path, [right, left]) == {p + 4: 2.0}

    def test_paired_end_fragment_contributes_both_ends(self, tmp_path):
        # Fragment [40, 140): R1 forward at 40, R2 reverse ending at 140.
        r1 = self._read("frag", 40, mate=110, read1=True)
        r2 = self._read("frag", 110, reverse=True, mate=40, read1=False)
        assert self._cuts(tmp_path, [r1, r2]) == {44: 1.0, 135: 1.0}

    def test_single_end_read_contributes_only_its_start(self, tmp_path):
        forward = self._read("f", 40)
        reverse = self._read("r", 110, reverse=True)
        # A reverse read "starts" at its 5' end: the right end of the alignment,
        # not the left-most coordinate.
        assert self._cuts(tmp_path, [forward]) == {44: 1.0}
        assert self._cuts(tmp_path, [reverse]) == {110 + self.READ - 5: 1.0}

    def test_base_quality_does_not_matter_but_filters_do(self, tmp_path):
        reads = [
            self._read("lowq", 40, baseq=2),  # counted: cuts ignore base quality
            self._read("dup", 60, flags=0x400),  # duplicate: dropped
            self._read("mapq", 80, mapq=5),  # below --min_mapq: dropped
        ]
        assert self._cuts(tmp_path, reads, min_mapq=20, min_baseq=30) == {44: 1.0}

    def test_cut_outside_the_region_is_not_counted(self, tmp_path):
        # Read overlaps [50, 100) but its cut (44) lies before it.
        assert self._cuts(tmp_path, [self._read("f", 40)], start=50, end=100) == {}

    def test_extend_size_broadcasts_each_cut(self, tmp_path):
        cuts = self._cuts(tmp_path, [self._read("f", 40)], extend_size=2)
        assert cuts == {42: 1.0, 43: 1.0, 44: 1.0, 45: 1.0, 46: 1.0}

    def test_run_bam2bw_writes_a_cut_track(self, tmp_path):
        bam = self._bam(
            tmp_path,
            [
                self._read("frag", 40, mate=110, read1=True),
                self._read("frag", 110, reverse=True, mate=40, read1=False),
            ],
        )
        run_bam2bw(
            bam_path=bam,
            fasta_path=self._fasta(tmp_path),
            out_dir=str(tmp_path / "o"),
            out_name="tn5",
            min_mapq=0,
            event="tn5",
        )
        with pyBigWig.open(str(tmp_path / "o" / "tn5.bw")) as bw:
            values = np.nan_to_num(np.array(bw.values("chr1", 0, self.LENGTH)))
        assert {int(i) for i in np.nonzero(values)[0]} == {44, 135}

    def test_ratio_mode_is_rejected(self, tmp_path):
        bam = self._bam(tmp_path, [self._read("f", 40)])
        with pytest.raises(ValueError, match="ratio"):
            run_bam2bw(
                bam_path=bam,
                fasta_path=self._fasta(tmp_path),
                out_dir=str(tmp_path),
                out_name="x",
                mode="ratio",
                event="tn5",
            )

    def test_unknown_event_is_rejected(self, tmp_path):
        bam = self._bam(tmp_path, [self._read("f", 40)])
        with pytest.raises(ValueError, match="event"):
            run_bam2bw(
                bam_path=bam,
                fasta_path=self._fasta(tmp_path),
                out_dir=str(tmp_path),
                out_name="x",
                event="cuts",
            )
