"""Build a deamination-aware BWA index for a reference FASTA.

Following the bwa-meth strategy, the index is built on a *doubly converted*
copy of the reference: every chromosome appears twice — once as a C->T
converted forward strand (prefixed ``f``) and once as a G->A converted
reverse-orientation strand (prefixed ``r``). This lets BWA-MEM map both
top-strand-derived (C->T pattern) and bottom-strand-derived (G->A pattern
in read orientation) deaminated reads against a single index.

Outputs (all next to the input FASTA):

  - ``<fasta>.fai``                                   — samtools faidx of the original
  - ``<fasta>.deamtools.c2t``                         — doubly-converted reference
  - ``<fasta>.deamtools.c2t.{amb,ann,bwt,pac,sa}``    — BWA-MEM index files
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


def run_index(fasta_path: str, force: bool = False) -> None:
    """Build the deamtools BWA index for ``fasta_path``.

    Skips work that has already been done unless ``force=True``.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    _check_executable("bwa")
    _check_executable("samtools")

    logger.info(f"Indexing {fasta_path}")

    fai_path = fasta_path + ".fai"
    if force or not os.path.exists(fai_path):
        logger.info("samtools faidx ...")
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
    else:
        logger.info(f"  {fai_path} exists; skipping faidx")

    converted_path = fasta_path + ".deamtools.c2t"
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

    logger.info("Done")
