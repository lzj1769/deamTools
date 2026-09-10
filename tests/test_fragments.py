"""Tests for deamtools.utils.fragments (mate pairing and merging)."""

import pysam
import pytest

from deamtools.utils import iter_fragments, mates_can_pair, merge_fragment_bases

REF_LEN = 20


def _read(
    name: str,
    seq: str,
    pos: int = 0,
    *,
    paired: bool = True,
    read1: bool = True,
    mate_unmapped: bool = False,
    mate_ref_id: int | None = 0,
    baseq: int = 40,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    flag = 0
    if paired:
        flag |= 0x1 | (0x40 if read1 else 0x80)
    if mate_unmapped:
        flag |= 0x8
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 30
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * len(seq))
    if mate_ref_id is not None:
        a.next_reference_id = mate_ref_id
        a.next_reference_start = pos
    return a


class TestMatesCanPair:
    def test_paired_read_with_mapped_mate_on_same_contig(self):
        assert mates_can_pair(_read("a", "ACGT")) is True

    def test_unpaired_read(self):
        assert mates_can_pair(_read("a", "ACGT", paired=False)) is False

    def test_mate_unmapped(self):
        assert mates_can_pair(_read("a", "ACGT", mate_unmapped=True)) is False

    def test_mate_on_another_contig(self):
        assert mates_can_pair(_read("a", "ACGT", mate_ref_id=1)) is False


class TestIterFragments:
    def _names(self, fragments):
        return [tuple(r.query_name for r in f) for f in fragments]

    def test_mates_are_paired_by_name(self):
        reads = [_read("a", "ACGT"), _read("a", "ACGT", read1=False)]
        assert self._names(iter_fragments(reads)) == [("a", "a")]

    def test_interleaved_pairs_are_matched_correctly(self):
        reads = [
            _read("a", "ACGT"),
            _read("b", "ACGT"),
            _read("a", "ACGT", read1=False),
            _read("b", "ACGT", read1=False),
        ]
        assert self._names(iter_fragments(reads)) == [("a", "a"), ("b", "b")]

    def test_unpaired_reads_pass_straight_through(self):
        reads = [_read(n, "ACGT", paired=False) for n in "abc"]
        assert self._names(iter_fragments(reads)) == [("a",), ("b",), ("c",)]

    def test_orphan_is_yielded_at_the_end_not_dropped(self):
        # 'b' never gets a partner -- e.g. its mate was filtered upstream.
        reads = [
            _read("a", "ACGT"),
            _read("b", "ACGT"),
            _read("a", "ACGT", read1=False),
        ]
        assert self._names(iter_fragments(reads)) == [("a", "a"), ("b",)]

    def test_every_input_record_comes_back_exactly_once(self):
        reads = [
            _read("a", "ACGT"),
            _read("a", "ACGT", read1=False),
            _read("b", "ACGT"),
            _read("c", "ACGT", paired=False),
        ]
        out = [r for frag in iter_fragments(reads) for r in frag]
        assert sorted(r.query_name for r in out) == ["a", "a", "b", "c"]

    def test_empty_input(self):
        assert list(iter_fragments([])) == []


class TestMergeFragmentBases:
    def test_single_read_maps_positions_to_bases(self):
        bases = merge_fragment_bases([_read("a", "ACGT", pos=5)], 0)
        assert bases == {5: "A", 6: "C", 7: "G", 8: "T"}

    def test_overlapping_mates_collapse_to_one_entry_per_position(self):
        frag = [_read("a", "ACGT", pos=0), _read("a", "ACGT", pos=0, read1=False)]
        assert merge_fragment_bases(frag, 0) == {0: "A", 1: "C", 2: "G", 3: "T"}

    def test_disjoint_mates_are_unioned(self):
        frag = [_read("a", "AC", pos=0), _read("a", "GT", pos=5, read1=False)]
        assert merge_fragment_bases(frag, 0) == {0: "A", 1: "C", 5: "G", 6: "T"}

    def test_higher_base_quality_wins_a_disagreement(self):
        loud = _read("a", "TTTT", pos=0, baseq=40)
        quiet = _read("a", "AAAA", pos=0, read1=False, baseq=2)
        assert merge_fragment_bases([loud, quiet], 0)[0] == "T"
        assert merge_fragment_bases([quiet, loud], 0)[0] == "T"

    def test_equal_quality_disagreement_is_dropped_as_ambiguous(self):
        """Neither mate is more credible, so the position is not called.

        Quality scores are binned on modern instruments (Q2/Q12/Q23/Q37), so
        equal-quality disagreements are common; resolving them by input order
        would just mean "read 1 always wins".
        """
        a = _read("a", "TTTT", pos=0, baseq=37)
        b = _read("a", "AAAA", pos=0, read1=False, baseq=37)
        assert merge_fragment_bases([a, b], 0) == {}
        assert merge_fragment_bases([b, a], 0) == {}

    def test_equal_quality_agreement_is_kept(self):
        a = _read("a", "ACGT", pos=0, baseq=37)
        b = _read("a", "ACGT", pos=0, read1=False, baseq=37)
        assert merge_fragment_bases([a, b], 0) == {0: "A", 1: "C", 2: "G", 3: "T"}

    def test_only_the_disagreeing_position_is_dropped(self):
        a = _read("a", "ACGT", pos=0, baseq=37)
        b = _read("a", "ATGT", pos=0, read1=False, baseq=37)
        assert merge_fragment_bases([a, b], 0) == {0: "A", 2: "G", 3: "T"}

    def test_min_baseq_drops_positions(self):
        frag = [_read("a", "ACGT", pos=0, baseq=5)]
        assert merge_fragment_bases(frag, 20) == {}

    def test_min_baseq_can_leave_only_the_confident_mate(self):
        loud = _read("a", "TTTT", pos=0, baseq=40)
        quiet = _read("a", "AAAA", pos=0, read1=False, baseq=5)
        assert merge_fragment_bases([loud, quiet], 20) == dict.fromkeys(range(4), "T")

    def test_region_bounds_are_applied(self):
        frag = [_read("a", "ACGTACGT", pos=0)]
        assert merge_fragment_bases(frag, 0, start=2, end=5) == {2: "G", 3: "T", 4: "A"}

    def test_read_without_a_sequence_is_skipped(self):
        empty = _read("a", "ACGT", pos=0)
        empty.query_sequence = None
        assert merge_fragment_bases([empty], 0) == {}

    @pytest.mark.parametrize("baseq", [0, 20, 40])
    def test_threshold_boundary_is_inclusive(self, baseq):
        frag = [_read("a", "ACGT", pos=0, baseq=baseq)]
        assert len(merge_fragment_bases(frag, baseq)) == 4
