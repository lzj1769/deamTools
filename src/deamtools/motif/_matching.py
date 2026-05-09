from __future__ import annotations
from collections.abc import Iterable
import MOODS.scan
import MOODS.tools

import argparse

from deamtools.utils import logger

def _get_motifs_from_jaspar(
    release: str = "JASPAR2024",
    collection: str = "CORE",
    tax_group: list[str] | None = None,
    all_versions: bool = False,
) -> Iterable | None:
    """
    Fetch transcription factor motifs from the JASPAR database.

    This function retrieves transcription factor binding motifs from the JASPAR database using the `pyjaspar` library.
    It allows filtering by JASPAR release, motif collection, taxonomic group, and version.

    Parameters
    ----------
    release :
        The release version ( e.g. JASPAR2020, JASPAR2024) of the JASPAR database to query.
    collection :
        The collection of motifs to query. Common options include:

            - `"CORE"`: High-quality, manually curated collection.
            - `"UNVALIDATED"`: Computationally predicted motifs.

    tax_group :
        A list of taxonomic groups to filter motifs. For example:

        - `["vertebrates"]`
        - `["plants", "insects"]`

        If `None`, defaults to `["vertebrates"]`.

    all_versions :
        Whether to fetch all versions of each motif. If `True`, retrieves all motif versions;
        otherwise, retrieves only the latest version.

    Returns
    -------
        An iterable of motif objects fetched from the JASPAR database.
        Returns `None` if the `pyjaspar` library is not installed or an error occurs.

    Raises
    ------
    ImportError
        If the `pyjaspar` library is not installed.

    Notes
    -----
    - Requires the `pyjaspar` library to interact with the JASPAR database. Install it via `pip install pyjaspar`.
    - The motifs fetched are represented as objects with attributes like `name`, `matrix_id`, and `counts`, which can be used for downstream analysis.

    Examples
    --------
    Fetch all motifs from the JASPAR2024 CORE collection for vertebrates:

    >>> motifs = get_motifs_from_jaspar(
    ...     release="JASPAR2024",
    ...     collection="CORE",
    ...     tax_group=["vertebrates"],
    ... )
    >>> print(len(motifs))
    """
    # check if JASPAR is installed
    try:
        from pyjaspar import jaspardb
    except ImportError:
        logger.error(
            "pyjaspar is not installed. Please install it first: pip install pyjaspar"
        )
        return None

    if tax_group is None:
        tax_group = ["vertebrates"]

    jdb_obj = jaspardb(release=release)
    motifs = jdb_obj.fetch_motifs(
        collection=collection, tax_group=tax_group, all_versions=all_versions
    )

    logger.info(f"Number of motifs fetched: {len(motifs)}")

    return motifs

def prepare_scanner(
    motifs: list,
    pseudocounts: float = 0.0001,
    p_value: float = 5e-05,
) -> MOODS.scan.Scanner:
    """
    Prepares a MOODS scanner for motif scanning.

    This function initializes a `MOODS.scan.Scanner` object with the given motifs,
    applying log-odds transformation to the motif counts and computing motif-specific
    score thresholds based on a given p-value. It also generates reverse complement
    motifs for scanning both DNA strands.


    Parameters
    ----------
    motifs :
         A list of motif objects containing nucleotide frequency counts
        (expected to have `.counts` attribute with keys "A", "C", "G", and "T").
    pseudocounts :
        A small pseudocount added to avoid zero probabilities in log-odds calculation.
    p_value :
        The statistical significance threshold for computing motif scoring thresholds.

    Returns
    -------
        A `MOODS.scan.Scanner` object initialized with the provided motifs and
        scoring thresholds, ready for scanning DNA sequences.

    Notes
    -----
        - The function creates both **original** and **reverse complement** motifs for scanning.
        - The background nucleotide frequency is assumed to be uniform (flat background).
        - `MOODS.tools.threshold_from_p()` is used to compute the score thresholds.

    Examples
    --------
    >>> import cell2net as cn
    >>> motifs = [some_motif_object1, some_motif_object2]  # Assume motif objects are loaded
    >>> scanner = cn.pp.prepare_scaner(motifs)
    >>> type(scanner)
    <class 'MOODS.scan.Scanner'>
    """
    n_motifs = len(motifs)
    bg = MOODS.tools.flat_bg(4)

    matrices = [None] * 2 * n_motifs
    thresholds = [None] * 2 * n_motifs
    for i, motif in enumerate(motifs):
        counts = (
            tuple(motif.counts["A"]),
            tuple(motif.counts["C"]),
            tuple(motif.counts["G"]),
            tuple(motif.counts["T"]),
        )

        matrices[i] = MOODS.tools.log_odds(counts, bg, pseudocounts)
        matrices[i + n_motifs] = MOODS.tools.reverse_complement(matrices[i])

        thresholds[i] = MOODS.tools.threshold_from_p(matrices[i], bg, p_value)
        thresholds[i + n_motifs] = thresholds[i]

    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices=matrices, bg=bg, thresholds=thresholds)

    return scanner

def run_motif_matching(genome_fasta: str = None,
                       bed_file: str = None,
                       output_file: str = None,
                       p_value_threshold: float = 1e-4) -> None:

    logger.info("Running motif matching...")
    logger.info(f"Genome FASTA: {genome_fasta}")
    logger.info(f"Input file: {bed_file}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"P-value threshold: {p_value_threshold}")



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