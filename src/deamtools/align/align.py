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

**The two mates are not converted differently.** Each candidate is *one*
conversion applied to the whole fragment in **reference** space; the labels
differ only because read 2's FASTQ sequence is the reverse complement of the
reference-space sequence, and the conversion has to be written in the space the
FASTQ is in. Complementing turns ``G->A`` into ``C->T``, so for any sequence
``s``::

    revcomp(s.replace("G", "A")) == revcomp(s).replace("C", "T")

which is exactly why read 2 carries the opposite label. That flip is what keeps
both mates on the *same* converted contig -- ``f`` sends both to ``f<chrom>``,
``r`` sends both to ``r<chrom>`` -- so bwa can pair them.

**Do not add the other two combinations** (both mates ``ct``, or both ``ga``).
They are not alternative hypotheses about the molecule: a ``ct`` mate maps to
``f<chrom>`` and a ``ga`` mate to ``r<chrom>``, so those combinations split the
pair across contigs. Measured over 300 simulated fragments, they put the mates
on one contig 7% of the time and produce **zero** proper pairs, yet still score
a mean ``AS(R1)+AS(R2)`` of ~177 -- above the valid ``r`` candidate's 159 -- so
feeding them to the take-best would sometimes select a non-pair over a correct
pair. See ``docs/algorithm.md`` for the full table.

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

import logging
import os
import re
import shutil
import subprocess
import threading
from collections.abc import Iterable
from typing import IO

import pysam

logger = logging.getLogger(__name__)

CT_TABLE = str.maketrans("Cc", "Tt")
GA_TABLE = str.maketrans("Gg", "Aa")
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# BWA-MEM options, matching bwa-meth 0.2.9 (`bwameth.py:bwa_mem`). A three-letter
# alignment leaves residual mismatches wherever the conversion runs the wrong way
# — for a dsDNA deaminase that is the minority-direction editing every read
# carries — so the defaults are miscalibrated for this data:
#   -T 40   raise the minimum output score (with -B 2 a good alignment scores high)
#   -B 2    halve the mismatch penalty; residual edits should not sink a read
#   -L 10   raise the clipping penalty, so bwa resolves rather than clips them
#   -M      mark shorter split hits as secondary
#   -C      carry the YS/YC comment through into the SAM (deamtools relies on it)
#   -U 100  (paired) heavily penalise unpaired placement
# Keeping these identical to bwa-meth's also makes the two directly comparable:
# any difference is then the read-conversion strategy, not parameter tuning.
BWA_MEM_OPTS = ("-T", "40", "-B", "2", "-L", "10", "-C", "-M")
BWA_MEM_PAIRED_OPTS = ("-U", "100")


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
    out: IO[str], name: str, seq: str, qual: str, candidate: str, original: str
) -> None:
    """Write one converted FASTQ record.

    The ``YS`` (original SEQ) and ``YC`` (candidate) tags are emitted as a
    tab-separated comment so ``bwa mem -C`` copies each as its own SAM tag.
    """
    out.write(f"@{name}\tYS:Z:{original}\tYC:Z:{candidate}\n{seq}\n+\n{qual}\n")


def _record_fields(record: pysam.FastxRecord) -> tuple[str, str, str]:
    """``(name, sequence, quality)`` of a FASTQ record, with types narrowed.

    pysam types all three as optional, because a FASTA record carries no
    quality and a truncated record may carry nothing. Anything fed to the
    aligner must have a name and a sequence, so a record missing either is
    reported here rather than surfacing as a ``NoneType`` error deep inside the
    conversion. A missing quality string is filled with ``I`` (Phred 40), which
    is what the caller did before.
    """
    name, seq = record.name, record.sequence
    if name is None or seq is None:
        raise ValueError(
            f"malformed FASTQ record: name={name!r}, sequence missing"
            if seq is None
            else f"malformed FASTQ record {name!r}: no name"
        )
    qual = record.quality if record.quality is not None else "I" * len(seq)
    return name, seq, qual


def _feed_converted(
    read1: str,
    read2: str | None,
    out: IO[str],
) -> None:
    """Stream both converted candidates of every read/fragment to ``out``."""
    try:
        if read2 is None:
            with pysam.FastxFile(read1) as fq:
                for r in fq:
                    name, seq, qual = _record_fields(r)
                    _write_record(out, name, seq.translate(CT_TABLE), qual, "ct", seq)
                    _write_record(out, name, seq.translate(GA_TABLE), qual, "ga", seq)
        else:
            with pysam.FastxFile(read1) as fq1, pysam.FastxFile(read2) as fq2:
                for pair_number, (r1, r2) in enumerate(
                    zip(fq1, fq2, strict=True), start=1
                ):
                    n1, s1, q1 = _record_fields(r1)
                    n2, s2, q2 = _record_fields(r2)
                    key1 = n1[:-2] if n1.endswith("/1") else n1
                    key2 = n2[:-2] if n2.endswith("/2") else n2
                    if key1 != key2:
                        raise ValueError(
                            f"FASTQ mate names differ at pair {pair_number}: "
                            f"{n1!r} != {n2!r}"
                        )
                    n1 = n2 = key1
                    # Orientation f: read1 C->T, read2 G->A (interleaved pair).
                    _write_record(out, n1, s1.translate(CT_TABLE), q1, "f", s1)
                    _write_record(out, n2, s2.translate(GA_TABLE), q2, "f", s2)
                    # Orientation r: read1 G->A, read2 C->T.
                    _write_record(out, n1, s1.translate(GA_TABLE), q1, "r", s1)
                    _write_record(out, n2, s2.translate(CT_TABLE), q2, "r", s2)
    finally:
        out.close()


def _emit_clean_header(fasta_path: str, out: IO[str]) -> None:
    """Write a fresh @HD + @SQ block from the original FASTA's .fai."""
    out.write("@HD\tVN:1.6\tSO:coordinate\n")
    with open(fasta_path + ".fai") as f:
        for line in f:
            chrom, length = line.split("\t")[:2]
            out.write(f"@SQ\tSN:{chrom}\tLN:{length}\n")


def _original_contig(name: str) -> str:
    return name[1:] if name.startswith(("f", "r")) else name


def _reference_tags(
    seq: str, cigar: str, chrom: str, pos: int, reference: pysam.FastaFile
) -> tuple[int, str, int, int]:
    """Return NM, MD, C→T and G→A counts against the original reference.

    Position is zero-based. Counts cover aligned bases only, without a base
    quality filter; sequence is already in reference orientation.
    """
    ops = [(int(n), op) for n, op in CIGAR_OP_RE.findall(cigar)]
    span = sum(n for n, op in ops if op in "MDN=X")
    ref = reference.fetch(chrom, pos, pos + span).upper()
    if len(ref) != span:
        raise ValueError(f"Alignment extends beyond reference: {chrom}:{pos}")
    seq = seq.upper()
    q = r = nm = matches = ct = ga = 0
    md: list[str] = []
    for n, op in ops:
        if op in "M=X":
            for i in range(n):
                if seq[q + i] == ref[r + i] and seq[q + i] in "ACGT":
                    matches += 1
                else:
                    md.extend((str(matches), ref[r + i]))
                    matches = 0
                    nm += 1
                    ct += ref[r + i] == "C" and seq[q + i] == "T"
                    ga += ref[r + i] == "G" and seq[q + i] == "A"
            q += n
            r += n
        elif op == "D":
            md.extend((str(matches), "^" + ref[r : r + n]))
            matches = 0
            nm += n
            r += n
        elif op == "N":
            r += n
        elif op == "I":
            nm += n
            q += n
        elif op == "S":
            q += n
    md.append(str(matches))
    return nm, "".join(md), ct, ga


def _restore_hit_tag(tag: str, original: str, reference: pysam.FastaFile) -> str:
    """Restore contigs and edit distances inside BWA SA/XA hit lists."""
    hits = []
    for hit in tag[5:].split(";"):
        if not hit:
            continue
        fields = hit.split(",")
        fields[0] = _original_contig(fields[0])
        if tag.startswith("SA:"):
            pos, strand, cigar = int(fields[1]) - 1, fields[2], fields[3]
        else:
            pos, strand, cigar = abs(int(fields[1])) - 1, fields[1][0], fields[2]
        seq = _revcomp(original) if strand == "-" else original
        left, right = _hard_clip_offsets(cigar)
        seq = seq[left : len(seq) - right]
        nm, _, _, _ = _reference_tags(seq, cigar, fields[0], pos, reference)
        fields[-1] = str(nm)
        hits.append(",".join(fields))
    return tag[:5] + ";".join(hits) + ";"


def _restore_alignment(line: str, reference: pysam.FastaFile | None = None) -> str:
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
    hit_tags: list[str] = []
    for tag in fields[11:]:
        if tag.startswith("YS:Z:"):
            orig = tag[5:]
        elif tag.startswith("YC:Z:"):
            continue  # candidate marker; internal only
        elif tag.startswith(("NM:", "MD:", "ZC:", "ZG:")):
            continue  # Converted-reference tags are invalid after SEQ restoration.
        elif tag.startswith(("SA:Z:", "XA:Z:")):
            hit_tags.append(tag)
        else:
            kept_tags.append(tag)

    original = orig
    if orig is not None and seq_field != "*":
        if flag & 16:
            orig = _revcomp(orig)
        left, right = _hard_clip_offsets(cigar)
        if left or right:
            orig = orig[left : len(orig) - right]
        if len(orig) != len(seq_field):
            raise ValueError(
                f"Cannot restore sequence for {fields[0]!r}: length mismatch"
            )
        fields[9] = orig

    if reference is not None and original is not None:
        kept_tags.extend(_restore_hit_tag(tag, original, reference) for tag in hit_tags)
        if not flag & 4 and seq_field != "*":
            nm, md, ct, ga = _reference_tags(
                fields[9], cigar, fields[2], int(fields[3]) - 1, reference
            )
            kept_tags.extend((f"NM:i:{nm}", f"MD:Z:{md}", f"ZC:i:{ct}", f"ZG:i:{ga}"))

    fields = fields[:11] + kept_tags
    return "\t".join(fields) + "\n"


def _tag_value(fields: list[str], prefix: str) -> str | None:
    """Value of a SAM tag (e.g. ``"YC:Z:"`` or ``"AS:i:"``) in ``fields[11:]``."""
    for tag in fields[11:]:
        if tag.startswith(prefix):
            return tag[len(prefix) :]
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


def _primary_locations(lines: list[str]) -> tuple[tuple, ...]:
    """Compare primary placements in original reference space, by mate."""
    locations = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & (0x100 | 0x800 | 0x4):
            continue
        locations.append((flag & 0xC0, _original_contig(f[2]), f[3], flag & 16, f[5]))
    return tuple(sorted(locations))


def _flush_group(
    lines: list[str], out: IO[str], reference: pysam.FastaFile | None = None
) -> None:
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
    best_locations = _primary_locations(by_candidate[best])
    competitors = [
        _primary_score(records)
        for key, records in by_candidate.items()
        if key != best
        and _primary_locations(records)
        and _primary_locations(records) != best_locations
    ]
    # Conservative score-gap ceiling, not a calibrated error probability.
    # Equal-scoring distinct placements get MAPQ 0; same-location conversion
    # duplicates do not reduce confidence. Never increase BWA's MAPQ.
    cap = (
        max(0, _primary_score(by_candidate[best]) - max(competitors))
        if competitors
        else None
    )
    for line in by_candidate[best]:
        fields = line.rstrip("\n").split("\t")
        if cap is not None:
            fields[4] = str(min(int(fields[4]), cap))
            # SA MAPQs describe the same selected candidate's split records.
            for i in range(11, len(fields)):
                if fields[i].startswith("SA:Z:"):
                    hits = []
                    for hit in fields[i][5:].rstrip(";").split(";"):
                        parts = hit.split(",")
                        parts[4] = str(min(int(parts[4]), cap))
                        hits.append(",".join(parts))
                    fields[i] = "SA:Z:" + ";".join(hits) + ";"
        out.write(_restore_alignment("\t".join(fields) + "\n", reference))


def _process_sam(
    bwa_stdout: Iterable[str],
    sort_stdin: IO[str],
    reference: pysam.FastaFile | None = None,
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
            _flush_group(group, sort_stdin, reference)
            group = []
        group_qname = qname
        group.append(line)
    if group:
        _flush_group(group, sort_stdin, reference)


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

    bwa_cmd = ["bwa", "mem", *BWA_MEM_OPTS, "-t", str(threads)]
    if paired:
        bwa_cmd += [*BWA_MEM_PAIRED_OPTS, "-p"]
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
        with open(sam_path, "w") as sam, pysam.FastaFile(fasta_path) as reference:
            _emit_clean_header(fasta_path, sam)
            if bwa_proc.stdout is None:  # unreachable with stdout=PIPE
                raise RuntimeError("bwa mem produced no stdout stream")
            _process_sam(bwa_proc.stdout, sam, reference)
    except BaseException:
        bwa_proc.kill()
        bwa_proc.wait()
        raise
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
