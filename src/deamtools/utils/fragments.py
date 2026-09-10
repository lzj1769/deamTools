"""Grouping paired-end records into fragments and merging their mates.

Whenever the insert is shorter than twice the read length the two mates overlap,
and every reference position in that overlap is reported twice -- once by each
mate. It is one position on one molecule, so counting it twice inflates coverage
and, because the overlap is the middle of the fragment rather than a random
subset of positions, biases any rate computed from it.

Both :mod:`deamtools.qc.qc` and :mod:`deamtools.preprocessing.bam2bw` therefore
tally per *fragment*: group records with :func:`iter_fragments`, collapse each
fragment's mates with :func:`merge_fragment_bases`, then count. Keeping that in
one place is deliberate -- the two modules previously drifted on exactly this
point. :mod:`deamtools.preprocessing.bam2fragment` predates these helpers and
merges mates itself, strand-aware, into a set of editing positions.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence

import pysam

# Stand-in quality for a record with no QUAL string at all ("*"). Such a record
# should not be silently dropped by a --min_baseq threshold, so it is treated as
# the SAM "quality unavailable" value rather than as quality 0.
_NO_QUAL = 255

# A fragment is the one or two records that came off the same molecule.
Fragment = tuple[pysam.AlignedSegment, ...]


def mates_can_pair(read: pysam.AlignedSegment) -> bool:
    """Could this record's mate turn up in the same single-contig pass?

    False for unpaired reads, for a read whose mate is unmapped, and for a read
    whose mate is on another contig -- in each case the record is a fragment on
    its own as far as this pass is concerned.
    """
    return (
        read.is_paired
        and not read.mate_is_unmapped
        and read.next_reference_id == read.reference_id
    )


def iter_fragments(reads: Iterable[pysam.AlignedSegment]) -> Iterator[Fragment]:
    """Group records into fragments, pairing mates by read name.

    Yields a 1- or 2-tuple per fragment. A record that could have a mate is held
    until its partner arrives; whatever is still waiting when ``reads`` runs out
    is yielded on its own, so **no input record is ever dropped** -- a mate can
    fail to arrive because it was filtered out upstream, not only because the
    file is odd.

    ``reads`` must already be filtered: everything it yields becomes part of a
    fragment. Memory is bounded by how many mates are in flight at once, which
    for a coordinate-sorted file is a function of the insert-size distribution;
    pass one contig at a time to bound it properly.
    """
    pending: dict[str, pysam.AlignedSegment] = {}
    for read in reads:
        if not mates_can_pair(read):
            yield (read,)
            continue
        qname = read.query_name or ""
        mate = pending.pop(qname, None)
        if mate is None:
            pending[qname] = read
            continue
        yield (mate, read)
    for orphan in pending.values():
        yield (orphan,)


def merge_fragment_bases(
    fragment: Sequence[pysam.AlignedSegment],
    min_baseq: int,
    start: int | None = None,
    end: int | None = None,
) -> dict[int, str]:
    """Collapse a fragment's mates into ``reference position -> called base``.

    Where the mates overlap they are reporting the same duplex: library prep
    turns a deaminated C into a real T:A pair, so both mates carry the event and
    are expected to agree. A disagreement there is therefore a sequencing error
    in one of them, and the higher base quality wins.

    Positions failing ``min_baseq`` are dropped, as are positions outside
    ``[start, end)`` when those are given. Only matched (M/=/X) bases contribute,
    so indels are skipped the way ``get_aligned_pairs(matches_only=True)``
    skips them.
    """
    best: dict[int, tuple[str, int]] = {}
    for read in fragment:
        seq = read.query_sequence
        if seq is None:
            continue
        quals = read.query_qualities
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if start is not None and ref_pos < start:
                continue
            if end is not None and ref_pos >= end:
                continue
            qual = quals[query_pos] if quals is not None else _NO_QUAL
            if qual < min_baseq:
                continue
            previous = best.get(ref_pos)
            if previous is None or qual > previous[1]:
                best[ref_pos] = (seq[query_pos], qual)
    return {ref_pos: base for ref_pos, (base, _) in best.items()}
