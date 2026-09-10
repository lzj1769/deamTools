"""Convert a coordinate-sorted BAM to a per-fragment editing-signal table.

Each output row represents a unique fragment defined by ``(chrom, start, end,
c_to_t_positions, g_to_a_positions[, barcode])``. The ``count`` field reports
the number of reads/pairs that produced that exact signature, and the two
editing columns are ``|``-separated lists of reference positions (0-based).

Output formats:

  * Without ``--barcode``::

        chrom  start  end  count  c2t_positions  g2a_positions

  * With ``--barcode`` (10x fragments-style ordering)::

        chrom  start  end  barcode  count  c2t_positions  g2a_positions

Either editing column is ``.`` when that fragment shows no event of that kind.

**The two edit directions are reported separately, and both are called
regardless of read orientation.** A reference ``C`` read as ``T`` means the top
strand was deaminated at that position; a reference ``G`` read as ``A`` means
the bottom strand was. Which strand carries the event is a property of the
molecule, not of the record that happens to observe it: a BAM stores every read
in reference orientation, and library prep fixes the deaminated U into a real
T:A pair, so both mismatch patterns are visible on reads of either orientation.

Until 2026-09-10 this module recorded ``C->T`` only on forward reads and
``G->A`` only on reverse reads, which used ``read.is_reverse`` as a proxy for
"which strand was edited". On ``data/ACCESS-ATAC/chr10.bam`` that proxy is only
about 65% accurate and the rule discarded **35.5%** of editing events -- events
carrying the same ``TC`` enzyme-motif fingerprint as the ones it kept, so
genuine deamination rather than noise. Keeping the two directions in separate
columns preserves the strand information the old rule was reaching for; the
ACCESS-ATAC preprint calls resolving both strands a route to better allelic
occupancy imputation.
"""

from __future__ import annotations

import gzip
import logging
import os
from collections import defaultdict
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import IO

import pysam

from deamtools.utils import get_chrom_sizes_from_bam, merge_fragment_bases

logger = logging.getLogger(__name__)

# Per-chromosome key layout: (start, end, c_to_t, g_to_a, barcode_or_None).
_FragKey = tuple[int, int, tuple[int, ...], tuple[int, ...], "str | None"]


def _open_output(path: str) -> IO[str]:
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def _editing_positions(
    fragment: Sequence[pysam.AlignedSegment],
    ref_seq: str,
    min_baseq: int,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Reference positions where this fragment records a deamination event.

    Returns ``(c_to_t, g_to_a)``: sorted positions where a reference ``C`` reads
    ``T`` (top-strand event) and where a reference ``G`` reads ``A``
    (bottom-strand event). Both are called **regardless of read orientation** --
    see the module docstring for why ``read.is_reverse`` is the wrong test.

    The mates are collapsed by :func:`~deamtools.utils.merge_fragment_bases`
    first, so a position both of them cover is one observation and the higher
    base quality settles a disagreement, rather than a low-quality mate being
    able to inject an event the confident mate contradicts.
    """
    bases = merge_fragment_bases(fragment, min_baseq)
    c_to_t: list[int] = []
    g_to_a: list[int] = []
    for ref_pos in sorted(bases):
        ref_base = ref_seq[ref_pos]
        read_base = bases[ref_pos]
        if ref_base == "C" and read_base == "T":
            c_to_t.append(ref_pos)
        elif ref_base == "G" and read_base == "A":
            g_to_a.append(ref_pos)
    return tuple(c_to_t), tuple(g_to_a)


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
    buffer: dict[str, pysam.AlignedSegment] = {}

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
                c_to_t, g_to_a = _editing_positions((read,), ref_seq, min_baseq)
                bc = _get_barcode(read, barcode_tag) if barcode else None
                counter[(start, end, c_to_t, g_to_a, bc)] += 1
                continue

            # Paired-end: only proper pairs whose mate is also mapped.
            if not read.is_proper_pair or read.mate_is_unmapped:
                continue

            qname = read.query_name
            if qname is None:  # unnamed record; cannot be paired up
                continue
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
            # Checked one by one rather than with `None in (...)`, which reads
            # the same but tells a type checker nothing about the four names.
            if r1_start is None or r2_start is None or r1_end is None or r2_end is None:
                continue

            start = min(r1_start, r2_start)
            end = max(r1_end, r2_end)

            c_to_t, g_to_a = _editing_positions((r1, r2), ref_seq, min_baseq)

            bc = None
            if barcode:
                bc = _get_barcode(r1, barcode_tag) or _get_barcode(r2, barcode_tag)

            counter[(start, end, c_to_t, g_to_a, bc)] += 1

    return chrom, dict(counter)


def _format_row(
    chrom: str,
    key: _FragKey,
    count: int,
    include_barcode: bool,
) -> str:
    start, end, c_to_t, g_to_a, barcode = key
    c2t_str = "|".join(str(p) for p in c_to_t) if c_to_t else "."
    g2a_str = "|".join(str(p) for p in g_to_a) if g_to_a else "."
    if include_barcode:
        bc_str = barcode if barcode is not None else "."
        return f"{chrom}\t{start}\t{end}\t{bc_str}\t{count}\t{c2t_str}\t{g2a_str}"
    return f"{chrom}\t{start}\t{end}\t{count}\t{c2t_str}\t{g2a_str}"


def run_bam2fragment(
    bam_path: str,
    fasta_path: str,
    out_dir: str,
    out_name: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
    threads: int = 1,
    barcode: bool = False,
    barcode_tag: str = "CB",
    gzip: bool = False,
) -> None:
    """Convert ``bam_path`` to a fragment table with per-fragment editing signals.

    The table is written to ``<out_dir>/<out_name>.tsv`` (or ``.tsv.gz`` when
    ``gzip=True``).
    """
    logger.info("Running bam2fragment")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")
    if barcode:
        logger.info(f"Barcode tag: {barcode_tag}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_sizes = get_chrom_sizes_from_bam(bam)
    chroms = list(chrom_sizes.keys())
    logger.info(f"Processing {len(chroms)} chromosome(s) with {threads} thread(s)")

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, f"{out_name}.tsv" + (".gz" if gzip else ""))

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
