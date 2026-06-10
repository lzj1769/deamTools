"""Build a deamination-aware BWA index for a reference FASTA.

Following the bwa-meth strategy, the index is built on a *doubly converted*
copy of the reference: every chromosome appears twice — once with all
cytosines converted to thymine (``C->T``, prefixed ``f``) and once with all
guanines converted to adenine (``G->A``, prefixed ``r``). Both are conversions
of the *forward* sequence; the ``f``/``r`` prefixes denote which deaminated
read population maps there (top-strand-derived ``C->T`` reads to ``f``,
bottom-strand-derived ``G->A`` reads to ``r``). Reducing the alphabet this way
lets BWA-MEM map heavily deaminated reads, and both mates of a pair land on the
same converted contig so proper pairing is preserved. The ``f``/``r`` prefix is
stripped from the chromosome name during ``deamtools align``.

Outputs:

  - ``<fasta>.fai``                                   — samtools faidx of the original
    (always written next to the FASTA; required by the pysam-based subcommands)
  - ``<out_dir>/<out_name>.deamtools.c2t``            — doubly-converted reference
  - ``<out_dir>/<out_name>.deamtools.c2t.{amb,ann,bwt,pac,sa}`` — BWA-MEM index files

``out_dir`` / ``out_name`` default to the FASTA's own directory and file name.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess

logger = logging.getLogger(__name__)

CT_TABLE = str.maketrans("Cc", "Tt")
GA_TABLE = str.maketrans("Gg", "Aa")
LINE_WIDTH = 80


def _check_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"'{name}' was not found in PATH. Install it before running 'deamtools index'."
        )


def _convert_fasta(fasta_path: str, output_path: str) -> None:
    """Stream FASTA, writing both ``f`` (C->T) and ``r`` (G->A) entries per chromosome."""

    def flush(header: str | None, seq_parts: list[str], out) -> None:
        if header is None:
            return
        seq = "".join(seq_parts)
        out.write(f">f{header}\n")
        for i in range(0, len(seq), LINE_WIDTH):
            out.write(seq[i : i + LINE_WIDTH].translate(CT_TABLE) + "\n")
        out.write(f">r{header}\n")
        for i in range(0, len(seq), LINE_WIDTH):
            out.write(seq[i : i + LINE_WIDTH].translate(GA_TABLE) + "\n")

    header: str | None = None
    seq_parts: list[str] = []

    with open(fasta_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                flush(header, seq_parts, fout)
                header = line[1:].split()[0].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush(header, seq_parts, fout)


def run_index(
    fasta_path: str,
    out_dir: str | None = None,
    out_name: str | None = None,
    force: bool = False,
) -> None:
    """Build the deamtools BWA index for ``fasta_path``.

    The converted reference and its BWA-MEM index are written to
    ``<out_dir>/<out_name>.deamtools.c2t*``. When ``out_dir`` / ``out_name`` are
    omitted they default to the FASTA's own directory and file name, so the
    index lands next to the FASTA — the location :func:`deamtools.align.run_align`
    looks in by default.

    The standard FASTA index (``<fasta>.fai``) is always written next to the
    original FASTA regardless of ``out_dir`` / ``out_name``, because the
    pysam-based subcommands (``bam2bw``, ``bam2fragment``, ``qc``, ...) require
    it there.

    Skips work that has already been done unless ``force=True``.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA to index.
    out_dir : str, optional
        Directory for the converted reference + BWA index. Defaults to the
        FASTA's directory.
    out_name : str, optional
        Base name for the converted reference + BWA index. Defaults to the
        FASTA file name.
    force : bool, default False
        Rebuild outputs even if they already exist.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    _check_executable("bwa")
    _check_executable("samtools")

    if out_dir is None:
        out_dir = os.path.dirname(fasta_path) or "."
    if out_name is None:
        out_name = os.path.basename(fasta_path)
    os.makedirs(out_dir, exist_ok=True)

    logger.info(f"Indexing {fasta_path}")

    # The .fai must sit next to the original FASTA for pysam-based subcommands.
    fai_path = fasta_path + ".fai"
    if force or not os.path.exists(fai_path):
        logger.info("samtools faidx ...")
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
    else:
        logger.info(f"  {fai_path} exists; skipping faidx")

    converted_path = os.path.join(out_dir, f"{out_name}.deamtools.c2t")
    if force or not os.path.exists(converted_path):
        logger.info(f"  writing converted reference to {converted_path}")
        _convert_fasta(fasta_path, converted_path)
    else:
        logger.info(f"  {converted_path} exists; skipping conversion")

    bwt_path = converted_path + ".bwt"
    if force or not os.path.exists(bwt_path):
        logger.info("bwa index (this may take a while) ...")
        subprocess.run(["bwa", "index", converted_path], check=True)
    else:
        logger.info(f"  {bwt_path} exists; skipping bwa index")

    logger.info(f"Converted index: {converted_path}")
    logger.info("Done")
