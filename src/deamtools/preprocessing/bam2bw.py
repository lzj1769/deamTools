from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pyBigWig
import pysam

from deamtools.utils import get_chrom_sizes_from_bam, get_chrom_sizes_from_file

logger = logging.getLogger(__name__)


def _load_regions(bed_path: str) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            regions.setdefault(chrom, []).append((start, end))
    # Merge overlapping intervals to avoid double-counting reads
    for chrom in regions:
        ivs = sorted(regions[chrom])
        merged: list[tuple[int, int]] = [ivs[0]]
        for start, end in ivs[1:]:
            if start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        regions[chrom] = merged
    return regions


def _count_deamination_on_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    chrom_size: int,
    regions: list[tuple[int, int]] | None,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
) -> tuple[str, np.ndarray]:
    counts = np.zeros(chrom_size, dtype=np.float32)

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        fetch_regions = regions if regions is not None else [(0, chrom_size)]

        for region_start, region_end in fetch_regions:
            region_end = min(region_end, chrom_size)
            ref_seq = fasta.fetch(chrom, region_start, region_end).upper()

            for read in bam.fetch(chrom, region_start, region_end):
                if (
                    read.is_unmapped
                    or read.is_duplicate
                    or read.is_qcfail
                    or read.is_secondary
                    or read.is_supplementary
                ):
                    continue
                if read.mapping_quality < min_mapq:
                    continue

                seq = read.query_sequence
                if seq is None:
                    continue

                quals = read.query_qualities
                is_reverse = read.is_reverse

                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if ref_pos < region_start or ref_pos >= region_end:
                        continue
                    if quals is not None and quals[query_pos] < min_baseq:
                        continue

                    ref_base = ref_seq[ref_pos - region_start]
                    read_base = seq[query_pos]

                    # Forward strand: C→T deamination
                    # Reverse strand: G→A (deamination of C on the template strand)
                    if is_reverse:
                        if ref_base == "G" and read_base == "A":
                            counts[ref_pos] += 1
                    else:
                        if ref_base == "C" and read_base == "T":
                            counts[ref_pos] += 1

    if extend_size > 0:
        kernel = np.ones(2 * extend_size + 1, dtype=np.float32)
        counts = np.convolve(counts, kernel, mode="same")

    return chrom, counts


def run_bam2bw(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    chrom_sizes_path: str | None = None,
    bed_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    extend_size: int = 0,
    threads: int = 1,
) -> None:
    logger.info("Running bam2bw")
    logger.info(f"BAM:   {bam_path}")
    logger.info(f"FASTA: {fasta_path}")

    if chrom_sizes_path is not None:
        chrom_sizes = get_chrom_sizes_from_file(chrom_sizes_path)
    else:
        logger.info("Inferring chromosome sizes from BAM header")
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            chrom_sizes = get_chrom_sizes_from_bam(bam)

    bed_regions: dict[str, list[tuple[int, int]]] | None = None
    if bed_path is not None:
        logger.info(f"Regions: {bed_path}")
        bed_regions = _load_regions(bed_path)
        logger.info(f"  {sum(len(v) for v in bed_regions.values())} intervals on "
                    f"{len(bed_regions)} chromosome(s)")

    chroms_to_process = [c for c in chrom_sizes if bed_regions is None or c in bed_regions]
    logger.info(f"Processing {len(chroms_to_process)} chromosome(s) "
                f"with {threads} thread(s)")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    results: dict[str, np.ndarray] = {}

    def _process(chrom: str) -> tuple[str, np.ndarray]:
        regions = bed_regions.get(chrom) if bed_regions is not None else None
        return _count_deamination_on_chrom(
            bam_path=bam_path,
            fasta_path=fasta_path,
            chrom=chrom,
            chrom_size=chrom_sizes[chrom],
            regions=regions,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            extend_size=extend_size,
        )

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_process, chrom): chrom for chrom in chroms_to_process}
        for future in as_completed(futures):
            chrom, counts = future.result()
            logger.info(f"  {chrom}: {int(counts.sum())} deamination event(s)")
            results[chrom] = counts

    logger.info(f"Writing {output_path}")
    with pyBigWig.open(output_path, "w") as bw:
        bw.addHeader(list(chrom_sizes.items()))
        for chrom in chroms_to_process:
            counts = results[chrom]
            nonzero = np.nonzero(counts)[0]
            if len(nonzero) == 0:
                continue
            bw.addEntries(
                chrom,
                nonzero.tolist(),
                values=counts[nonzero].tolist(),
                span=1,
            )

    logger.info("Done")
