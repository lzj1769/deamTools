"""Tests for deamtools.preprocessing.bam2bw."""

import os

import pyBigWig
import pysam
import pytest

from deamtools.preprocessing.bam2bw import (
    _count_deamination_on_chrom,
    _load_regions,
    run_bam2bw,
)

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


class TestLoadRegions:
    def test_basic(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\n")
        assert _load_regions(str(bed)) == {"chr1": [(0, 5)]}

    def test_overlapping_intervals_merged(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\nchr1\t3\t8\n")
        assert _load_regions(str(bed)) == {"chr1": [(0, 8)]}

    def test_non_overlapping_intervals_kept(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t3\nchr1\t5\t8\n")
        assert _load_regions(str(bed)) == {"chr1": [(0, 3), (5, 8)]}

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("# header\nchr1\t0\t5\n")
        assert _load_regions(str(bed)) == {"chr1": [(0, 5)]}

    def test_multiple_chromosomes(self, tmp_path):
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t5\nchr2\t10\t20\n")
        result = _load_regions(str(bed))
        assert set(result.keys()) == {"chr1", "chr2"}


# ---------------------------------------------------------------------------
# _count_deamination_on_chrom
# ---------------------------------------------------------------------------


class TestCountDeamination:
    """Unit tests for the per-chromosome counting kernel."""

    def _run(self, bam_path, fasta_file, *, regions=None, min_mapq=0, min_baseq=0, extend_size=0):
        _, counts = _count_deamination_on_chrom(
            bam_path=bam_path,
            fasta_path=fasta_file,
            chrom="chr1",
            chrom_size=len(REF_SEQ),
            regions=regions,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
        )
        return counts

    # -- Deamination detection ------------------------------------------------

    def test_forward_ct_event_detected(self, tmp_path, fasta_file):
        # REF:  A C G T C G A T C G
        # READ: A T G T C G A T C G  ← C→T at pos 1
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        counts = self._run(bam, fasta_file)
        assert counts[1] == 1.0
        assert counts.sum() == 1.0

    def test_reverse_ga_event_detected(self, tmp_path, fasta_file):
        # Reverse-strand read; G→A at ref pos 5 (- strand C deaminated to T).
        # query_sequence derived by reverse-complementing the deaminated - strand:
        #   normal rev-comp of REF = "CGATCGACGT"
        #   with deamination at ref pos 5: raw read = "CGATTGACGT"
        #   stored in BAM (rev-comp of raw) = "ACGTCAATCG"
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ACGTCAATCG", 0, 0, is_reverse=True)])
        counts = self._run(bam, fasta_file)
        assert counts[5] == 1.0
        assert counts.sum() == 1.0

    def test_no_event_on_exact_match(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", REF_SEQ, 0, 0)])
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
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTTGATCG", 0, 0)])
        counts = self._run(bam, fasta_file)
        assert counts[1] == 1.0
        assert counts[4] == 1.0
        assert counts.sum() == 2.0

    # -- Filtering ------------------------------------------------------------

    def test_low_mapq_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, mapq=10)])
        assert self._run(bam, fasta_file, min_mapq=20).sum() == 0.0

    def test_mapq_at_threshold_included(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, mapq=20)])
        assert self._run(bam, fasta_file, min_mapq=20).sum() == 1.0

    def test_low_baseq_position_excluded(self, tmp_path, fasta_file):
        # baseq=5 < min_baseq=20 → the C→T event at pos 1 is not counted
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, baseq=5)])
        assert self._run(bam, fasta_file, min_baseq=20).sum() == 0.0

    def test_secondary_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x100)])
        assert self._run(bam, fasta_file).sum() == 0.0

    def test_duplicate_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x400)])
        assert self._run(bam, fasta_file).sum() == 0.0

    def test_supplementary_read_excluded(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0, extra_flags=0x800)])
        assert self._run(bam, fasta_file).sum() == 0.0

    # -- Region restriction ---------------------------------------------------

    def test_event_inside_region_counted(self, tmp_path, fasta_file):
        # C→T at pos 1 and pos 4; restrict to [0, 3) → only pos 1 counted
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTTGATCG", 0, 0)])
        counts = self._run(bam, fasta_file, regions=[(0, 3)])
        assert counts[1] == 1.0
        assert counts[4] == 0.0
        assert counts.sum() == 1.0

    def test_event_outside_region_excluded(self, tmp_path, fasta_file):
        # C→T at pos 1; restrict to [5, 10) → no events in region
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        counts = self._run(bam, fasta_file, regions=[(5, 10)])
        assert counts.sum() == 0.0

    # -- Signal extension -----------------------------------------------------

    def test_extend_size_spreads_signal(self, tmp_path, fasta_file):
        # Single C→T at pos 1; extend_size=2 → signal at positions 0,1,2,3
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        counts = self._run(bam, fasta_file, extend_size=2)
        assert counts[0] == 1.0
        assert counts[1] == 1.0
        assert counts[2] == 1.0
        assert counts[3] == 1.0
        assert counts[4] == 0.0

    def test_extend_size_zero_unchanged(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        counts = self._run(bam, fasta_file, extend_size=0)
        assert counts[1] == 1.0
        assert counts.sum() == 1.0


# ---------------------------------------------------------------------------
# run_bam2bw (integration)
# ---------------------------------------------------------------------------


class TestRunBam2bw:
    """End-to-end tests that write a real BigWig and check its content."""

    def test_bigwig_created_with_correct_signal(self, tmp_path, fasta_file, chrom_sizes_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        out = str(tmp_path / "out.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=chrom_sizes_file, min_mapq=0, min_baseq=0)

        assert os.path.exists(out)
        with pyBigWig.open(out) as bw:
            # Event at pos 1
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(1.0)
            # No event at pos 0
            assert bw.stats("chr1", 0, 1, type="mean")[0] is None

    def test_infer_chrom_sizes_from_bam(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTCGATCG", 0, 0)])
        out = str(tmp_path / "out.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=None, min_mapq=0, min_baseq=0)
        assert os.path.exists(out)

    def test_bed_restricts_output_signal(self, tmp_path, fasta_file, chrom_sizes_file):
        # C→T at pos 1 and pos 4; BED restricts to [0, 3)
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", "ATGTTGATCG", 0, 0)])
        bed = tmp_path / "r.bed"
        bed.write_text("chr1\t0\t3\n")
        out = str(tmp_path / "out.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=chrom_sizes_file, bed_path=str(bed),
                   min_mapq=0, min_baseq=0)

        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(1.0)
            assert bw.stats("chr1", 4, 5, type="mean")[0] is None

    def test_output_parent_directory_created(self, tmp_path, fasta_file, chrom_sizes_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", REF_SEQ, 0, 0)])
        out = str(tmp_path / "nested" / "dir" / "out.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=chrom_sizes_file, min_mapq=0, min_baseq=0)
        assert os.path.exists(out)

    def test_no_events_writes_empty_bigwig(self, tmp_path, fasta_file, chrom_sizes_file):
        # Read identical to reference → no deamination events
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", REF_SEQ, 0, 0)])
        out = str(tmp_path / "out.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=chrom_sizes_file, min_mapq=0, min_baseq=0)

        assert os.path.exists(out)
        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 0, len(REF_SEQ), type="mean")[0] is None


# ---------------------------------------------------------------------------
# Ratio mode
# ---------------------------------------------------------------------------


class TestRatioMode:
    """Tests for the --mode ratio output: events / informative coverage."""

    def _run(self, bam_path, fasta_file, *, mode, extend_size=0):
        _, signal = _count_deamination_on_chrom(
            bam_path=bam_path,
            fasta_path=fasta_file,
            chrom="chr1",
            chrom_size=len(REF_SEQ),
            regions=None,
            min_mapq=0,
            min_baseq=0,
            extend_size=extend_size,
            mode=mode,
        )
        return signal

    def test_invalid_mode_raises(self, tmp_path, fasta_file):
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", REF_SEQ, 0, 0)])
        with pytest.raises(ValueError):
            _count_deamination_on_chrom(
                bam_path=bam, fasta_path=fasta_file, chrom="chr1",
                chrom_size=len(REF_SEQ), regions=None,
                min_mapq=0, min_baseq=0, extend_size=0, mode="bogus",
            )

    def test_ratio_three_quarters_at_C(self, tmp_path, fasta_file):
        # Reference C at pos 1: 3 reads with T (events), 1 read with C (no event).
        # Expected ratio at pos 1 = 3 / (3 + 1) = 0.75.
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at pos 1 -> event
            _make_read("r2", "ATGTCGATCG", 0, 0),  # T at pos 1 -> event
            _make_read("r3", "ATGTCGATCG", 0, 0),  # T at pos 1 -> event
            _make_read("r4", REF_SEQ,      0, 0),  # C at pos 1 -> coverage only
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file, mode="ratio")
        assert ratio[1] == pytest.approx(0.75)
        # No events at any other C/G position -> ratio 0 (or coverage 0).
        # Position 4 (ref C) is covered by all 4 reads with C -> 0/4 = 0.
        assert ratio[4] == 0.0

    def test_ratio_uninformative_base_excluded_from_denominator(self, tmp_path, fasta_file):
        # At ref C (pos 1): 1 read T (event), 1 read A (uninformative).
        # Denominator only counts C/T reads -> 1 / (1 + 0) = 1.0, not 0.5.
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at pos 1
            _make_read("r2", "AAGTCGATCG", 0, 0),  # A at pos 1 (uninformative)
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file, mode="ratio")
        assert ratio[1] == pytest.approx(1.0)

    def test_ratio_zero_when_no_coverage(self, tmp_path, fasta_file):
        # Read identical to reference -> events 0 everywhere; positions with
        # coverage have ratio 0; uncovered positions also 0.
        bam = _write_bam(str(tmp_path / "t.bam"),
                         [_make_read("r1", REF_SEQ, 0, 0)])
        ratio = self._run(bam, fasta_file, mode="ratio")
        assert (ratio == 0.0).all()

    def test_ratio_reverse_strand(self, tmp_path, fasta_file):
        # Reverse-strand reads at ref G (pos 2). Read with A at pos 2 -> event.
        # Expected ratio: 2 events / 3 informative reads.
        seq = "ACATCGATCG"  # A at pos 2 (deam event under reverse strand)
        reads = [
            _make_read("r1", seq,     0, 0, is_reverse=True),  # A
            _make_read("r2", seq,     0, 0, is_reverse=True),  # A
            _make_read("r3", REF_SEQ, 0, 0, is_reverse=True),  # G (coverage only)
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file, mode="ratio")
        assert ratio[2] == pytest.approx(2.0 / 3.0)

    def test_ratio_with_extend_size_uses_window_sums(self, tmp_path, fasta_file):
        # Two reads: one event at pos 1 (T at C), one read with C at pos 1.
        # extend_size=1 -> events convolved over [0,1,2], coverage convolved
        # over [0,1,2]. At pos 0,1,2: ratio = (sum events in window) / (sum cov in window).
        # events: [0, 1, 0, ...]; coverage at C-positions only:
        #   pos 1: 2 (one T, one C); other positions: 0
        # After convolve with window=3:
        #   events: [1, 1, 1, 0, ...]
        #   coverage: [2, 2, 2, 0, ...]
        #   ratio at 0,1,2: 0.5 each.
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at C (event)
            _make_read("r2", REF_SEQ,      0, 0),  # C at C (coverage only)
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        ratio = self._run(bam, fasta_file, mode="ratio", extend_size=1)
        for i in (0, 1, 2):
            assert ratio[i] == pytest.approx(0.5)

    def test_ratio_bigwig_round_trip(self, tmp_path, fasta_file, chrom_sizes_file):
        reads = [
            _make_read("r1", "ATGTCGATCG", 0, 0),  # T at C (event)
            _make_read("r2", "ATGTCGATCG", 0, 0),  # T at C (event)
            _make_read("r3", REF_SEQ,      0, 0),  # C at C (coverage only)
        ]
        bam = _write_bam(str(tmp_path / "t.bam"), reads)
        out = str(tmp_path / "ratio.bw")
        run_bam2bw(bam_path=bam, fasta_path=fasta_file, output_path=out,
                   chrom_sizes_path=chrom_sizes_file,
                   min_mapq=0, min_baseq=0, mode="ratio")
        with pyBigWig.open(out) as bw:
            assert bw.stats("chr1", 1, 2, type="mean")[0] == pytest.approx(2.0 / 3.0)
