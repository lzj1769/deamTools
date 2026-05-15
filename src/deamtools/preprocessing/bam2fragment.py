"""Convert a coordinate-sorted BAM to a per-fragment editing-signal table.

Each output row represents a unique fragment defined by ``(chrom, start, end,
editing_positions[, barcode])``. The ``count`` field reports the number of
reads/pairs that produced that exact signature, and ``editing_positions`` is
a ``|``-separated list of reference positions (0-based) where the fragment
shows a C->T (forward read) or G->A (reverse read) deamination event.

Output formats:

  * Without ``--barcode``::

        chrom  start  end  count  pos1|pos2|...

  * With ``--barcode`` (10x fragments-style ordering)::

        chrom  start  end  barcode  count  pos1|pos2|...

Fragments with no detected editing positions emit ``.`` in the editing column.
"""

from __future__ import annotations

import gzip
import logging
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import IO

import pysam

from deamtools.utils import get_chrom_sizes_from_bam

logger = logging.getLogger(__name__)

# Per-chromosome key layout: (start, end, edits, barcode_or_None).
_FragKey = tuple[int, int, tuple[int, ...], "str | None"]


def _open_output(path: str) -> IO[str]:
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def _editing_positions(read, ref_seq: str, min_baseq: int) -> list[int]:
    """Reference positions (0-based) where this read records a deamination event."""
    seq = read.query_sequence
    if seq is None:
        return []
    quals = read.query_qualities
    is_reverse = read.is_reverse
    positions: list[int] = []
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if quals is not None and quals[query_pos] < min_baseq:
            continue
        ref_base = ref_seq[ref_pos]
        read_base = seq[query_pos]
        if is_reverse:
            if ref_base == "G" and read_base == "A":
                positions.append(ref_pos)
        else:
            if ref_base == "C" and read_base == "T":
                positions.append(ref_pos)
    return positions


def _passes_basic_filters(read) -> bool:
    return not (
        read.is_unmapped
        or read.is_duplicate
        or read.is_qcfail
        or read.is_secondary
        or read.is_supplementary
    )


def _get_barcode(read, tag: str) -> str | None:
    try:
        return read.get_tag(tag)
    except KeyError:
        return None


def _process_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    min_mapq: int,
    min_baseq: int,
    barcode: bool,
    barcode_tag: str,
) -> tuple[str, dict[_FragKey, int]]:
    counter: dict[_FragKey, int] = defaultdict(int)
    buffer: dict[str, "pysam.AlignedSegment"] = {}

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        ref_seq = fasta.fetch(chrom).upper()

        for read in bam.fetch(chrom):
            if not _passes_basic_filters(read):
                continue

            if not read.is_paired:
                if read.mapping_quality < min_mapq:
                    continue
                start = read.reference_start
                end = read.reference_end
                if start is None or end is None:
                    continue
                positions = tuple(_editing_positions(read, ref_seq, min_baseq))
                bc = _get_barcode(read, barcode_tag) if barcode else None
                counter[(start, end, positions, bc)] += 1
                continue

            # Paired-end: only proper pairs whose mate is also mapped.
            if not read.is_proper_pair or read.mate_is_unmapped:
                continue

            qname = read.query_name
            mate = buffer.pop(qname, None)
            if mate is None:
                buffer[qname] = read
                continue

            # Apply MAPQ at pair completion so a low-MAPQ mate drops the pair.
            if mate.mapping_quality < min_mapq or read.mapping_quality < min_mapq:
                continue

            r1, r2 = (mate, read) if mate.is_read1 else (read, mate)
            r1_start, r2_start = r1.reference_start, r2.reference_start
            r1_end, r2_end = r1.reference_end, r2.reference_end
            if None in (r1_start, r2_start, r1_end, r2_end):
                continue

            start = min(r1_start, r2_start)
            end = max(r1_end, r2_end)

            edits: set[int] = set()
            edits.update(_editing_positions(r1, ref_seq, min_baseq))
            edits.update(_editing_positions(r2, ref_seq, min_baseq))
            positions = tuple(sorted(edits))

            bc: str | None = None
            if barcode:
                bc = _get_barcode(r1, barcode_tag) or _get_barcode(r2, barcode_tag)

            counter[(start, end, positions, bc)] += 1

    return chrom, dict(counter)


def _format_row(
    chrom: str,
    key: _FragKey,
    count: int,
    include_barcode: bool,
) -> str:
    start, end, positions, barcode = key
    pos_str = "|".join(str(p) for p in positions) if positions else "."
    if include_barcode:
        bc_str = barcode if barcode is not None else "."
        return f"{chrom}\t{start}\t{end}\t{bc_str}\t{count}\t{pos_str}"
    return f"{chrom}\t{start}\t{end}\t{count}\t{pos_str}"


def run_bam2fragment(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
    threads: int = 1,
    barcode: bool = False,
    barcode_tag: str = "CB",
) -> None:
    """Convert ``bam_path`` to a fragment file with per-fragment editing signals."""
    logger.info("Running bam2fragment")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")
    if barcode:
        logger.info(f"Barcode tag: {barcode_tag}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_sizes = get_chrom_sizes_from_bam(bam)
    chroms = list(chrom_sizes.keys())
    logger.info(f"Processing {len(chroms)} chromosome(s) with {threads} thread(s)")

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    results: dict[str, dict[_FragKey, int]] = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(
                _process_chrom,
                bam_path,
                fasta_path,
                chrom,
                min_mapq,
                min_baseq,
                barcode,
                barcode_tag,
            ): chrom
            for chrom in chroms
        }
        for future in as_completed(futures):
            chrom, counter = future.result()
            logger.info(
                f"  {chrom}: {len(counter)} unique fragment signature(s) "
                f"({sum(counter.values())} total)"
            )
            results[chrom] = counter

    logger.info(f"Writing {output_path}")
    with _open_output(output_path) as out:
        for chrom in chroms:
            counter = results.get(chrom, {})
            for key in sorted(counter):
                out.write(_format_row(chrom, key, counter[key], barcode) + "\n")

    logger.info("Done")
