"""Motif matching against a reference genome with MOODS.

Scans the sequence of a set of genomic regions (e.g. accessible peaks) for
occurrences of transcription-factor motifs and writes the hits as a BED file of
motif-predicted binding sites (MPBSs). Scanning uses MOODS: log-odds matrices,
p-value-derived score thresholds, and reverse-complement matrices for both
strands.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterable

import MOODS.scan
import MOODS.tools
import pysam

from deamtools.utils import _load_regions

logger = logging.getLogger(__name__)


def _get_motifs_from_jaspar(
    release: str = "JASPAR2024",
    collection: str = "CORE",
    tax_group: list[str] | None = None,
    all_versions: bool = False,
) -> Iterable | None:
    """Fetch transcription-factor motifs from the JASPAR database.

    Retrieves motifs via the optional ``pyjaspar`` library, filtered by release,
    collection, and taxonomic group.

    Parameters
    ----------
    release : str
        JASPAR release (e.g. ``"JASPAR2020"``, ``"JASPAR2024"``).
    collection : str
        Motif collection, e.g. ``"CORE"`` (curated) or ``"UNVALIDATED"``.
    tax_group : list[str], optional
        Taxonomic groups to keep (default ``["vertebrates"]``).
    all_versions : bool
        Fetch every motif version instead of only the latest.

    Returns
    -------
    Iterable or None
        Motif objects (each with ``.counts``, ``.matrix_id``, ``.name``), or
        ``None`` if ``pyjaspar`` is not installed.
    """
    try:
        from pyjaspar import jaspardb
    except ImportError:
        logger.error(
            "pyjaspar is not installed. Install it first: pip install pyjaspar"
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
    """Build a MOODS scanner for a list of motifs.

    Each motif's count matrix is converted to a log-odds matrix against a flat
    background, a score threshold is derived from ``p_value``, and the reverse
    complement is added so both strands are scanned. The matrices are laid out
    as ``[fwd_0, ..., fwd_{n-1}, rc_0, ..., rc_{n-1}]``, so
    ``scanner.scan(seq)[i]`` holds the forward hits and ``[i + n]`` the
    reverse-complement hits for motif ``i``.

    Parameters
    ----------
    motifs : list
        Motif objects exposing ``.counts`` with keys ``"A"``, ``"C"``, ``"G"``,
        ``"T"`` (each a per-position count sequence).
    pseudocounts : float
        Pseudocount added in the log-odds transform.
    p_value : float
        Significance threshold used to derive per-motif score cutoffs.

    Returns
    -------
    MOODS.scan.Scanner
        A scanner ready for :func:`scan_sequence`.
    """
    n_motifs = len(motifs)
    bg = MOODS.tools.flat_bg(4)

    matrices = [None] * (2 * n_motifs)
    thresholds = [None] * (2 * n_motifs)
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


def _motif_width(motif) -> int:
    return len(motif.counts["A"])


def _motif_name(motif) -> str:
    """A BED-friendly label: ``<matrix_id>.<name>`` when both are present."""
    matrix_id = getattr(motif, "matrix_id", None)
    name = getattr(motif, "name", None)
    if matrix_id and name and name != matrix_id:
        return f"{matrix_id}.{name}"
    return str(matrix_id or name or "motif")


def scan_sequence(
    scanner: MOODS.scan.Scanner,
    motifs: list,
    seq: str,
    chrom: str,
    offset: int = 0,
) -> list[tuple[str, int, int, str, float, str]]:
    """Scan one sequence and return motif matches as BED-style tuples.

    Parameters
    ----------
    scanner : MOODS.scan.Scanner
        Scanner from :func:`prepare_scanner` for the same ``motifs`` list.
    motifs : list
        The motifs the scanner was built from (defines order and width).
    seq : str
        DNA sequence to scan (upper-case A/C/G/T).
    chrom : str
        Chromosome name for the emitted coordinates.
    offset : int
        Genomic coordinate of ``seq[0]`` (e.g. the region start), added to every
        match position.

    Returns
    -------
    list[tuple]
        ``(chrom, start, end, name, score, strand)`` per match, with half-open
        ``[start, end)`` genomic coordinates.
    """
    results = scanner.scan(seq)
    n_motifs = len(motifs)
    matches: list[tuple[str, int, int, str, float, str]] = []
    for i, motif in enumerate(motifs):
        width = _motif_width(motif)
        name = _motif_name(motif)
        for strand, hits in (("+", results[i]), ("-", results[i + n_motifs])):
            for hit in hits:
                start = offset + int(hit.pos)
                matches.append(
                    (chrom, start, start + width, name, float(hit.score), strand)
                )
    return matches


def run_motif_matching(
    fasta_path: str,
    bed_path: str,
    output_path: str,
    motifs: list | None = None,
    release: str = "JASPAR2024",
    collection: str = "CORE",
    tax_group: list[str] | None = None,
    pseudocounts: float = 0.0001,
    p_value: float = 1e-4,
) -> None:
    """Scan BED regions for motif matches and write a BED of binding sites.

    For each interval in ``bed_path`` the reference sequence is read from
    ``fasta_path`` and scanned with MOODS; every hit above the ``p_value``
    threshold is written to ``output_path`` as a 6-column BED line
    ``chrom  start  end  motif  score  strand``.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA indexed with ``samtools faidx``.
    bed_path : str
        BED file of regions to scan (overlapping intervals are merged).
    output_path : str
        Output BED path. Parent directories are created.
    motifs : list, optional
        Pre-loaded motif objects. When ``None``, motifs are fetched from JASPAR
        using ``release`` / ``collection`` / ``tax_group`` (needs ``pyjaspar``).
    release, collection, tax_group : str / list[str]
        JASPAR query parameters used when ``motifs`` is ``None``.
    pseudocounts : float
        Pseudocount for the log-odds transform.
    p_value : float
        Significance threshold for motif hits.
    """
    logger.info("Running motif matching")
    logger.info(f"FASTA:   {fasta_path}")
    logger.info(f"Regions: {bed_path}")

    if motifs is None:
        motifs = _get_motifs_from_jaspar(
            release=release, collection=collection, tax_group=tax_group
        )
        if not motifs:
            raise RuntimeError(
                "No motifs available. Install pyjaspar (pip install pyjaspar) "
                "or pass motifs= explicitly."
            )
    motifs = list(motifs)
    logger.info(f"Motifs:  {len(motifs)} (p-value {p_value})")

    scanner = prepare_scanner(motifs, pseudocounts=pseudocounts, p_value=p_value)
    regions = _load_regions(bed_path)
    logger.info(f"Scanning {len(regions)} region(s)")

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    n_matches = 0
    with pysam.FastaFile(fasta_path) as fasta, open(output_path, "w") as out:
        for chrom, start, end in zip(
            regions["chrom"], regions["start"], regions["end"], strict=True
        ):
            chrom, start, end = str(chrom), int(start), int(end)
            seq = fasta.fetch(chrom, start, end).upper()
            if not seq:
                continue
            for c, s, e, name, score, strand in scan_sequence(
                scanner, motifs, seq, chrom, start
            ):
                out.write(f"{c}\t{s}\t{e}\t{name}\t{score:.4f}\t{strand}\n")
                n_matches += 1

    logger.info(f"Wrote {n_matches} match(es) to {output_path}")
