"""Deamination-aware alignment with a dual-conversion, take-best read converter.

  FASTQ(s) --[convert x2]--> bwa mem -C --[group + pick best]--> <out_name>.sam
      --[samtools sort]--> coordinate-sorted <out_name>.bam (+ .bai)

Because the ACCESS-ATAC deaminase edits cytosines on **both** strands, a single
read can carry both ``C->T`` and ``G->A`` deamination events, so converting it
in only one direction leaves the other as mismatches. Instead, every read is
emitted in **two** converted forms and bwa maps both; the better-scoring one is
kept:

* Single-end — two candidates per read: ``C->T`` (``YC:Z:ct``) and ``G->A``
  (``YC:Z:ga``).
* Paired-end — two fragment orientations, with a consistent direction for the
  whole pair: ``f`` = (read1 ``C->T``, read2 ``G->A``) and ``r`` =
  (read1 ``G->A``, read2 ``C->T``).

Both candidates of a read/fragment share the original read name; the candidate
is marked with a ``YC:Z:`` tag and the original sequence stashed in ``YS:Z:``
(both carried through ``bwa mem -C``). In post-processing, records are grouped
by read name and the candidate with the higher primary alignment score (sum of
the mates' ``AS`` for pairs) is kept; the original SEQ is restored from ``YS``,
the ``f``/``r`` prefix is stripped from RNAME/RNEXT, and the ``YS``/``YC`` tags
are dropped. The restored SAM is written to ``<out_name>.sam`` and converted to
a coordinate-sorted, indexed BAM with samtools.
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


def _write_record(
    out: io.TextIOBase, name: str, seq: str, qual: str, candidate: str, original: str
) -> None:
    """Write one converted FASTQ record.

    The ``YS`` (original SEQ) and ``YC`` (candidate) tags are emitted as a
    tab-separated comment so ``bwa mem -C`` copies each as its own SAM tag.
    """
    out.write(f"@{name}\tYS:Z:{original}\tYC:Z:{candidate}\n{seq}\n+\n{qual}\n")


def _feed_converted(
    read1: str,
    read2: str | None,
    out: io.TextIOBase,
) -> None:
    """Stream both converted candidates of every read/fragment to ``out``."""
    try:
        if read2 is None:
            with pysam.FastxFile(read1) as fq:
                for r in fq:
                    seq = r.sequence
                    qual = r.quality if r.quality is not None else "I" * len(seq)
                    _write_record(out, r.name, seq.translate(CT_TABLE), qual, "ct", seq)
                    _write_record(out, r.name, seq.translate(GA_TABLE), qual, "ga", seq)
        else:
            with pysam.FastxFile(read1) as fq1, pysam.FastxFile(read2) as fq2:
                for r1, r2 in zip(fq1, fq2, strict=True):
                    s1, s2 = r1.sequence, r2.sequence
                    q1 = r1.quality if r1.quality is not None else "I" * len(s1)
                    q2 = r2.quality if r2.quality is not None else "I" * len(s2)
                    # Orientation f: read1 C->T, read2 G->A (interleaved pair).
                    _write_record(out, r1.name, s1.translate(CT_TABLE), q1, "f", s1)
                    _write_record(out, r2.name, s2.translate(GA_TABLE), q2, "f", s2)
                    # Orientation r: read1 G->A, read2 C->T.
                    _write_record(out, r1.name, s1.translate(GA_TABLE), q1, "r", s1)
                    _write_record(out, r2.name, s2.translate(CT_TABLE), q2, "r", s2)
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
        elif tag.startswith("YC:Z:"):
            continue  # candidate marker; internal only
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


def _tag_value(fields: list[str], prefix: str) -> str | None:
    """Value of a SAM tag (e.g. ``"YC:Z:"`` or ``"AS:i:"``) in ``fields[11:]``."""
    for tag in fields[11:]:
        if tag.startswith(prefix):
            return tag[len(prefix):]
    return None


def _primary_score(lines: list[str]) -> int:
    """Sum of ``AS`` over a candidate's primary records (mates of a pair).

    Secondary (0x100) and supplementary (0x800) records are ignored; an
    unmapped primary contributes -1 so any mapped placement outranks it.
    """
    total = 0
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        if flag & 0x100 or flag & 0x800:
            continue
        if flag & 0x4:
            total += -1
            continue
        as_val = _tag_value(fields, "AS:i:")
        total += int(as_val) if as_val is not None else 0
    return total


def _flush_group(lines: list[str], out: io.TextIOBase) -> None:
    """Pick the best candidate among ``lines`` (one read name) and emit it.

    Records are partitioned by their ``YC`` candidate tag; the candidate with
    the highest primary alignment score is restored and written, the others are
    dropped. On a tie the first-seen candidate wins (``ct``/``f``).
    """
    by_candidate: dict[str, list[str]] = {}
    for line in lines:
        key = _tag_value(line.rstrip("\n").split("\t"), "YC:Z:") or ""
        by_candidate.setdefault(key, []).append(line)

    best = max(by_candidate, key=lambda k: _primary_score(by_candidate[k]))
    for line in by_candidate[best]:
        out.write(_restore_alignment(line))


def _process_sam(
    bwa_stdout: io.TextIOBase,
    sort_stdin: io.TextIOBase,
) -> None:
    """Group bwa output by read name and write the best candidate of each.

    bwa-mem preserves input order, so a read's two converted candidates (which
    share the read name) are emitted consecutively; records are buffered until
    the read name changes, then the best candidate is chosen and restored.
    """
    group: list[str] = []
    group_qname: str | None = None
    for line in bwa_stdout:
        if line.startswith("@"):
            # @HD and @SQ are emitted from the FASTA index; pass through
            # everything else (@PG, @RG, @CO).
            if line.startswith(("@HD", "@SQ")):
                continue
            sort_stdin.write(line)
            continue
        tab = line.find("\t")
        qname = line[:tab] if tab != -1 else line.rstrip("\n")
        if group and qname != group_qname:
            _flush_group(group, sort_stdin)
            group = []
        group_qname = qname
        group.append(line)
    if group:
        _flush_group(group, sort_stdin)


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
    sam_path = os.path.join(out_dir, f"{out_name}.sam")
    paired = read2 is not None
    logger.info(f"Aligning {'paired' if paired else 'single'}-end reads")
    logger.info(f"  R1:    {read1}")
    if paired:
        logger.info(f"  R2:    {read2}")
    logger.info(f"  Index: {converted_path}")
    logger.info(f"  Out:   {output_bam}")

    os.makedirs(out_dir, exist_ok=True)

    bwa_cmd = ["bwa", "mem", "-C", "-t", str(threads)]
    if paired:
        bwa_cmd += ["-p"]
    if read_group is not None:
        bwa_cmd += ["-R", read_group]
    bwa_cmd += [converted_path, "-"]

    # Step 1: bwa mem -> restore original sequences/names -> <out_name>.sam.
    logger.info(f"bwa mem ({threads}t), dual conversion + take-best -> {sam_path}")
    bwa_proc = subprocess.Popen(
        bwa_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
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
        with open(sam_path, "w") as sam:
            _emit_clean_header(fasta_path, sam)
            _process_sam(bwa_proc.stdout, sam)
    finally:
        feeder.join()

    bwa_rc = bwa_proc.wait()
    if feeder_exc:
        raise feeder_exc[0]
    if bwa_rc != 0:
        raise subprocess.CalledProcessError(bwa_rc, bwa_cmd)

    # Step 2: convert the SAM to a coordinate-sorted, indexed BAM.
    logger.info(f"samtools sort ({threads}t) {sam_path} -> {output_bam}")
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", output_bam, sam_path],
        check=True,
    )
    logger.info("samtools index ...")
    subprocess.run(["samtools", "index", output_bam], check=True)
    logger.info("Done")
