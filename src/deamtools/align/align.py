"""Deamination-aware alignment with a bwa-meth-style read converter.

The pipeline is:

  FASTQ(s) --[convert]--> bwa mem -C --[post-process]--> samtools sort --> BAM

Read 1 is converted C->T and read 2 (if paired) is converted G->A on the fly,
with the original read sequence stashed in a ``YS:Z:`` tag carried through
``bwa mem -C``. After alignment the original SEQ is restored from ``YS`` and
the ``f``/``r`` prefix is stripped from RNAME/RNEXT to recover original
chromosome names.
"""

from __future__ import annotations

import io
import logging
import os
import re
import shutil
import subprocess
import threading

import pysam

logger = logging.getLogger(__name__)

CT_TABLE = str.maketrans("Cc", "Tt")
GA_TABLE = str.maketrans("Gg", "Aa")
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _check_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"'{name}' was not found in PATH.")


def _revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def _hard_clip_offsets(cigar: str) -> tuple[int, int]:
    if "H" not in cigar:
        return 0, 0
    ops = CIGAR_OP_RE.findall(cigar)
    left = int(ops[0][0]) if ops and ops[0][1] == "H" else 0
    right = int(ops[-1][0]) if ops and ops[-1][1] == "H" else 0
    return left, right


def _write_converted(out: io.TextIOBase, entry, table) -> None:
    seq = entry.sequence
    qual = entry.quality if entry.quality is not None else "I" * len(seq)
    converted = seq.translate(table)
    # Tab-separated comment so that `bwa mem -C` emits each token as its own SAM tag.
    out.write(f"@{entry.name}\tYS:Z:{seq}\n{converted}\n+\n{qual}\n")


def _feed_converted(
    read1: str,
    read2: str | None,
    out: io.TextIOBase,
) -> None:
    try:
        if read2 is None:
            with pysam.FastxFile(read1) as fq:
                for entry in fq:
                    _write_converted(out, entry, CT_TABLE)
        else:
            with pysam.FastxFile(read1) as fq1, pysam.FastxFile(read2) as fq2:
                for r1, r2 in zip(fq1, fq2, strict=True):
                    _write_converted(out, r1, CT_TABLE)
                    _write_converted(out, r2, GA_TABLE)
    finally:
        out.close()


def _emit_clean_header(fasta_path: str, out: io.TextIOBase) -> None:
    """Write a fresh @HD + @SQ block from the original FASTA's .fai."""
    out.write("@HD\tVN:1.6\tSO:coordinate\n")
    with open(fasta_path + ".fai") as f:
        for line in f:
            chrom, length = line.split("\t")[:2]
            out.write(f"@SQ\tSN:{chrom}\tLN:{length}\n")


def _restore_alignment(line: str) -> str:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return line

    rname = fields[2]
    if rname not in ("", "*") and rname[0] in ("f", "r"):
        fields[2] = rname[1:]
    rnext = fields[6]
    if rnext not in ("", "*", "=") and rnext[0] in ("f", "r"):
        fields[6] = rnext[1:]

    flag = int(fields[1])
    cigar = fields[5]
    seq_field = fields[9]

    orig: str | None = None
    kept_tags: list[str] = []
    for tag in fields[11:]:
        if tag.startswith("YS:Z:"):
            orig = tag[5:]
        else:
            kept_tags.append(tag)

    if orig is not None and seq_field != "*":
        if flag & 16:
            orig = _revcomp(orig)
        left, right = _hard_clip_offsets(cigar)
        if left or right:
            orig = orig[left : len(orig) - right]
        if len(orig) == len(seq_field):
            fields[9] = orig

    fields = fields[:11] + kept_tags
    return "\t".join(fields) + "\n"


def _process_sam(
    bwa_stdout: io.TextIOBase,
    sort_stdin: io.TextIOBase,
) -> None:
    for line in bwa_stdout:
        if line.startswith("@"):
            # @HD and @SQ are emitted from the FASTA index; pass through
            # everything else (@PG, @RG, @CO).
            if line.startswith(("@HD", "@SQ")):
                continue
            sort_stdin.write(line)
            continue
        sort_stdin.write(_restore_alignment(line))


def run_align(
    fasta_path: str,
    read1: str,
    out_dir: str,
    out_name: str,
    read2: str | None = None,
    threads: int = 1,
    read_group: str | None = None,
    index_path: str | None = None,
) -> None:
    """Align deaminated reads and write a sorted, indexed BAM.

    The BAM is written to ``<out_dir>/<out_name>.bam`` (with a companion
    ``.bai`` index). The reference must already have been prepared with
    :func:`deamtools.align.index.run_index`.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA previously indexed with ``deamtools index``.
    read1 : str
        FASTQ for read 1 (or the only FASTQ for single-end input). Plain or
        gzipped.
    out_dir : str
        Output directory. Created if it does not exist.
    out_name : str
        Base name (without extension) for the output; writes
        ``<out_dir>/<out_name>.bam``.
    read2 : str, optional
        FASTQ for read 2 (paired-end). Omit for single-end alignment.
    threads : int, default 1
        Total threads, split between ``bwa mem`` and ``samtools sort``.
    read_group : str, optional
        Read-group line passed to ``bwa mem -R``.
    index_path : str, optional
        Path to the converted reference built by ``deamtools index``
        (``<out_dir>/<out_name>.deamtools.c2t``). Use this when the index was
        built with a custom ``--out_dir`` / ``--out_name``. Defaults to
        ``<fasta>.deamtools.c2t`` (next to the FASTA).
    """
    converted_path = (
        index_path if index_path is not None else fasta_path + ".deamtools.c2t"
    )
    if not os.path.exists(converted_path + ".bwt"):
        raise FileNotFoundError(
            f"BWA index not found at {converted_path}.bwt — run "
            f"'deamtools index --fasta {fasta_path}' first, and pass its "
            f"--out_dir/--out_name location here via --index if it was custom."
        )
    if not os.path.exists(fasta_path + ".fai"):
        raise FileNotFoundError(
            f"FASTA index not found at {fasta_path}.fai — "
            f"run 'deamtools index --fasta {fasta_path}' first."
        )
    if not os.path.exists(read1):
        raise FileNotFoundError(f"FASTQ not found: {read1}")
    if read2 is not None and not os.path.exists(read2):
        raise FileNotFoundError(f"FASTQ not found: {read2}")
    _check_executable("bwa")
    _check_executable("samtools")

    output_bam = os.path.join(out_dir, f"{out_name}.bam")
    paired = read2 is not None
    logger.info(f"Aligning {'paired' if paired else 'single'}-end reads")
    logger.info(f"  R1:    {read1}")
    if paired:
        logger.info(f"  R2:    {read2}")
    logger.info(f"  Index: {converted_path}")
    logger.info(f"  Out:   {output_bam}")

    os.makedirs(out_dir, exist_ok=True)

    bwa_threads = max(1, (threads + 1) // 2)
    sort_threads = max(1, threads - bwa_threads)

    bwa_cmd = ["bwa", "mem", "-C", "-t", str(bwa_threads)]
    if paired:
        bwa_cmd += ["-p"]
    if read_group is not None:
        bwa_cmd += ["-R", read_group]
    bwa_cmd += [converted_path, "-"]
    sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-o", output_bam, "-"]

    logger.info(f"  bwa mem ({bwa_threads}t) | post-process | samtools sort ({sort_threads}t)")

    bwa_proc = subprocess.Popen(
        bwa_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )

    feeder_exc: list[BaseException] = []

    def _feed():
        try:
            _feed_converted(read1, read2, bwa_proc.stdin)
        except BaseException as e:
            feeder_exc.append(e)

    feeder = threading.Thread(target=_feed, daemon=True)
    feeder.start()

    try:
        _emit_clean_header(fasta_path, sort_proc.stdin)
        _process_sam(bwa_proc.stdout, sort_proc.stdin)
    finally:
        sort_proc.stdin.close()
    feeder.join()

    bwa_rc = bwa_proc.wait()
    sort_rc = sort_proc.wait()

    if feeder_exc:
        raise feeder_exc[0]
    if bwa_rc != 0:
        raise subprocess.CalledProcessError(bwa_rc, bwa_cmd)
    if sort_rc != 0:
        raise subprocess.CalledProcessError(sort_rc, sort_cmd)

    logger.info("samtools index ...")
    subprocess.run(["samtools", "index", output_bam], check=True)
    logger.info("Done")
