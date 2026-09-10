"""Motif matching against a reference genome with motifmatchpy.

Scans the sequence of a set of genomic regions (e.g. accessible peaks) for
occurrences of transcription-factor motifs and writes the hits as a BED file of
motif-predicted binding sites (MPBSs). Scanning uses `motifmatchpy
<https://github.com/lzj1769/motifmatchpy>`_: log-odds matrices, p-value-derived
score thresholds, and both strands.

motifmatchpy replaced MOODS here. Its scanner takes both strands itself and
returns hits already carrying a motif name, an end coordinate and a strand, so
the reverse-complement matrices, the ``[fwd_0..fwd_n, rc_0..rc_n]`` layout and
the index arithmetic that recovered a hit's motif and strand are all gone.
Scores and coordinates are unchanged: both use the same log-odds construction
and report positions on the forward strand.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterable, Sequence

import motifmatchpy as mm
import pysam

from deamtools.utils import _load_regions

logger = logging.getLogger(__name__)

# Column order expected by the log-odds transform.
_BASES = ("A", "C", "G", "T")


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


def _motif_name(motif) -> str:
    """A BED-friendly label: ``<matrix_id>.<name>`` when both are present."""
    matrix_id = getattr(motif, "matrix_id", None)
    name = getattr(motif, "name", None)
    if matrix_id and name and name != matrix_id:
        return f"{matrix_id}.{name}"
    return str(matrix_id or name or "motif")


def _to_motif(motif, bg: list[float], pseudocounts: float) -> mm.Motif:
    """Coerce a motif into a ``motifmatchpy.Motif``.

    An ``mm.Motif`` (e.g. from :func:`load_motifs_from_files`) already holds a
    log-odds matrix and is returned unchanged. Anything else is treated as a
    JASPAR-style count matrix: any object exposing ``.counts`` with ``"A"``/
    ``"C"``/``"G"``/``"T"`` keys, which is what ``pyjaspar`` returns, whose
    counts are turned into a log-odds matrix against ``bg``.
    """
    if isinstance(motif, mm.Motif):
        return motif
    counts = tuple(tuple(motif.counts[base]) for base in _BASES)
    # log_odds returns a list of lists; Motif's annotation asks for tuples.
    matrix = tuple(tuple(row) for row in mm.tools.log_odds(counts, bg, pseudocounts))
    return mm.Motif(_motif_name(motif), matrix)


def load_motifs_from_files(
    paths: Sequence[str], pseudocounts: float = 0.0001
) -> list[mm.Motif]:
    """Read motifs from PFM/ADM files.

    Each motif is named after its file stem, so ``MA0001.1.pfm`` becomes
    ``MA0001.1`` in the BED ``name`` column. motifmatchpy would otherwise use
    the full filename, extension included.

    Parameters
    ----------
    paths : sequence of str
        Motif files. ``.pfm`` is a position frequency matrix (order 0);
        ``.adm`` is an adjacent dinucleotide model (order 1).
    pseudocounts : float
        Pseudocount added in the log-odds transform, matching what
        :func:`prepare_scanner` applies to JASPAR motifs.

    Returns
    -------
    list[motifmatchpy.Motif]
        Motifs ready to pass to :func:`prepare_scanner` or
        :func:`run_motif_matching`.

    Raises
    ------
    FileNotFoundError
        If any path does not exist. Checked up front, so a typo in a long list
        fails before any parsing work.
    """
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(f"motif file(s) not found: {', '.join(missing)}")

    names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
    motifs = mm.read_motifs(
        list(paths),
        bg=mm.tools.flat_bg(4),
        pseudocount=pseudocounts,
        names=names,
    )
    logger.info(f"Loaded {len(motifs)} motif(s) from {len(paths)} file(s)")
    return motifs


def prepare_scanner(
    motifs: list,
    pseudocounts: float = 0.0001,
    p_value: float = 5e-05,
) -> mm.MotifScanner:
    """Build a motif scanner for a list of motifs.

    Each motif's count matrix is converted to a log-odds matrix against a flat
    background and a score threshold is derived from ``p_value``. The scanner
    searches both strands.

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
    motifmatchpy.MotifScanner
        A scanner ready for :func:`scan_sequence`.
    """
    bg = mm.tools.flat_bg(4)
    return mm.MotifScanner(
        [_to_motif(m, bg, pseudocounts) for m in motifs],
        p_value=p_value,
        bg=bg,
        both_strands=True,
    )


def scan_sequence(
    scanner: mm.MotifScanner,
    seq: str,
    chrom: str,
    offset: int = 0,
) -> list[tuple[str, int, int, str, float, str]]:
    """Scan one sequence and return motif matches as BED-style tuples.

    Parameters
    ----------
    scanner : motifmatchpy.MotifScanner
        Scanner from :func:`prepare_scanner`.
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
        ``[start, end)`` genomic coordinates. Positions are on the forward
        strand for both strands' hits, so a minus-strand hit spans the same
        interval its plus-strand counterpart would.
    """
    return [
        (
            chrom,
            offset + int(hit.pos),
            offset + int(hit.end),
            hit.name,
            float(hit.score),
            hit.strand,
        )
        for hit in scanner.scan(seq, sequence_name=chrom)
    ]


def run_motif_matching(
    fasta_path: str,
    bed_path: str,
    out_dir: str,
    out_name: str,
    motifs: Sequence | None = None,
    motif_files: Sequence[str] | None = None,
    release: str = "JASPAR2024",
    collection: str = "CORE",
    tax_group: list[str] | None = None,
    pseudocounts: float = 0.0001,
    p_value: float = 1e-4,
) -> None:
    """Scan BED regions for motif matches and write a BED of binding sites.

    For each interval in ``bed_path`` the reference sequence is read from
    ``fasta_path`` and scanned; every hit above the ``p_value`` threshold is
    written to ``<out_dir>/<out_name>.bed`` as a 6-column BED line
    ``chrom  start  end  motif  score  strand``.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA indexed with ``samtools faidx``.
    bed_path : str
        BED file of regions to scan (overlapping intervals are merged).
    out_dir : str
        Output directory. Created if it does not exist.
    out_name : str
        Base name (without extension) for the output; writes
        ``<out_dir>/<out_name>.bed``.
    motifs : list, optional
        Pre-loaded motifs, either ``motifmatchpy.Motif`` objects or JASPAR-style
        count matrices. Ignored when ``motif_files`` is given.
    motif_files : sequence of str, optional
        PFM/ADM files to read motifs from, as an alternative to querying
        JASPAR. Takes precedence over ``motifs``.
    release, collection, tax_group : str / list[str]
        JASPAR query parameters, used only when neither ``motif_files`` nor
        ``motifs`` is given.
    pseudocounts : float
        Pseudocount for the log-odds transform.
    p_value : float
        Significance threshold for motif hits.
    """
    logger.info("Running motif matching")
    logger.info(f"FASTA:   {fasta_path}")
    logger.info(f"Regions: {bed_path}")

    if motif_files:
        resolved: Iterable = load_motifs_from_files(
            motif_files, pseudocounts=pseudocounts
        )
    elif motifs is not None:
        resolved = motifs
    else:
        fetched = _get_motifs_from_jaspar(
            release=release, collection=collection, tax_group=tax_group
        )
        if not fetched:
            raise RuntimeError(
                "No motifs available. Pass --motif_files, or install pyjaspar "
                "(pip install pyjaspar) to fetch them from JASPAR."
            )
        resolved = fetched
    motifs = list(resolved)
    logger.info(f"Motifs:  {len(motifs)} (p-value {p_value})")

    scanner = prepare_scanner(motifs, pseudocounts=pseudocounts, p_value=p_value)
    regions = _load_regions(bed_path)
    logger.info(f"Scanning {len(regions)} region(s)")

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, f"{out_name}.bed")

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
                scanner, seq, chrom, start
            ):
                out.write(f"{c}\t{s}\t{e}\t{name}\t{score:.4f}\t{strand}\n")
                n_matches += 1

    logger.info(f"Wrote {n_matches} match(es) to {output_path}")
