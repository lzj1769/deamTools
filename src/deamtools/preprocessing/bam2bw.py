from __future__ import annotations

import argparse
import logging

import pyBigWig

from deamtools.utils import get_chrom_sizes_from_bam, get_chrom_sizes_from_file

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def run_bam2bw(bam_path: str = None, 
               fasta_path: str = None,
               bed_path: str = None,
               chrom_sizes: dict[str, int] = None,
               output_path: str = None) -> None:

    logging.info("Running bam2bw...")
    logging.info(f"Input BAM: {bam_path}")
    logging.info(f"Reference FASTA: {fasta_path}")
    
    

    if chrom_sizes is None:
        logging.info("No chromosome sizes provided, inferring from BAM header...")
        with pyBigWig.open(bam_path) as bam_file:
            chrom_sizes = get_chrom_sizes_from_bam(bam_file)
    else:
        logging.info("Using provided chromosome sizes.")
        chrom_sizes = get_chrom_sizes_from_file(chrom_sizes)

    logging.info(f"Filtering regions: {bed_path if bed_path else 'None (whole genome)'}")
    logging.info(f"Output BigWig: {output_path}")


    logging.info(f"Total of {len(grs)} regions")

    if bed_path is not None:
        logging.info(f"Loading regions from BED file: {bed_path}")
        chrom_sizes = {}
        with open(bed_path, "r") as bed_file:
            for line in bed_file:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), end)

    with pyBigWig.open(output_path, "wb") as bw:
        bw.addHeader(list(chrom_sizes.items()))


    raise NotImplementedError("bam2bw functionality is not yet implemented.")


    logging.info("Done!")