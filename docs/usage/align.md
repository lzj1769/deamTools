# align

Align deaminated sequencing reads (single- or paired-end) to a reference prepared with [`deamtools index`](index.md), writing a coordinate-sorted, indexed BAM.

## Synopsis

```
deamtools align --fasta FILE --fastq1 FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--fasta FILE` | Reference FASTA. Must already have been indexed with `deamtools index` (i.e. `<fasta>.fai` and `<fasta>.deamtools.c2t.*` exist). |
| `--fastq1 FILE` | FASTQ for read 1, or the only FASTQ for single-end input. Plain or gzipped. |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. Writes a sorted, indexed `<out_dir>/<out_name>.bam` (and `.bam.bai`). |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--fastq2 FILE` | *(single-end)* | FASTQ for read 2. Provide it for paired-end alignment; omit for single-end. |
| `--index FILE` | *(next to the FASTA)* | Path to the converted reference built by `deamtools index` (`<out_dir>/<out_name>.deamtools.c2t`). Use this when the index was built with a custom `--out_dir`/`--out_name`. |
| `--read_group STR` | *(none)* | Read-group line passed to `bwa mem -R`, e.g. `'@RG\tID:s1\tSM:sample1\tLB:lib1\tPL:ILLUMINA'`. |
| `--threads INT` | `1` | Total threads, split between `bwa mem` and `samtools sort`. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

:::{note}
By default `align` looks for the converted index next to the FASTA (`<fasta>.deamtools.c2t`). If you built it with a custom location — `deamtools index --fasta ref.fa --out_dir idx --out_name foo` — point `align` at it with `--index idx/foo.deamtools.c2t`. The `index` command logs this exact path (`Converted index: …`) so you can copy it.
:::

## Requirements

`bwa` and `samtools` must be on your `PATH`, and the reference must have been prepared with `deamtools index` first.

## How it works

```
FASTQ(s) ──[convert]──▶ bwa mem -C ──[restore]──▶ samtools sort ──▶ BAM (+ .bai)
   R1: C→T                 |  three-letter        |  strip f/r prefix
   R2: G→A                 |  alignment to the     |  restore original SEQ
   + stash original (YS)   |  f/r converted ref    |  rebuild header
```

The whole pipeline runs as a stream — reads are converted, mapped, and restored on the fly and piped straight into `samtools sort` — so no intermediate files are written.

### Why the reads are converted

A deaminated read carries many `C→T` (or `G→A`) substitutions relative to the genome. Aligned directly, these look like dense mismatches and the read either maps with a poor score or fails to map, biasing against the most accessible (most edited) loci. DeamTools removes the bias by collapsing the alphabet — the same **three-letter (bisulfite-style) alignment** used by bwa-meth: if both the reads and the reference are `C→T`-converted, a genuine deamination event is no longer a mismatch, so editing density no longer affects mappability.

### 1. On-the-fly read conversion

Each read is converted to the three-letter space as it is streamed to the aligner:

- **Read 1** is `C→T`-converted (every C and c → T/t).
- **Read 2** (paired-end only) is `G→A`-converted — read 2 is sequenced from the complementary strand, where top-strand C→T deamination appears as G→A.

The **original, unconverted sequence is preserved** in a `YS:Z:` SAM tag so it can be restored after mapping. This is done by appending the tag as a tab-separated comment on the FASTQ header line and running `bwa mem -C`, which copies the comment verbatim into the SAM record:

```
FASTQ record written to bwa:        Resulting (pre-restore) SAM fields:
@read1  YS:Z:ACGTTCGA               SEQ  = ATGTTTGA        (C→T converted)
ATGTTTGA                            tag  = YS:Z:ACGTTCGA   (original, via -C)
+
IIIIIIII                            (base qualities are passed through unchanged)
```

For paired-end input the two FASTQs are read in lockstep and written **interleaved** (R1, R2, R1, R2, …); `bwa mem -p` then pairs consecutive records. Reads with no quality string (e.g. FASTA input) get a placeholder quality.

### 2. Mapping to the doubly-converted reference

`bwa mem -C -t <n> [-p] [-R <rg>] <converted_ref> -` reads the interleaved stream from stdin and maps against the reference built by [`deamtools index`](index.md), which contains, for every chromosome:

- `f<chrom>` — the `C→T`-converted forward sequence, and
- `r<chrom>` — the `G→A`-converted forward sequence.

A `C→T`-converted read 1 matches the `f` contigs; a `G→A`-converted read 2 matches the `r` contigs. Crucially, **both mates of a fragment map to the same converted contig** (the `f`/`r` copy for the strand the fragment derives from), so BWA still flags them as a proper pair and computes insert sizes correctly.

### 3. Restoring the alignments

BWA's SAM stream is rewritten line by line back into the original reference space before it reaches `samtools sort`:

- **Header.** BWA's `@SQ`/`@HD` lines describe the doubled `f`/`r` contigs, so they are dropped and replaced with a fresh `@HD VN:1.6 SO:coordinate` plus one `@SQ` per chromosome read from the original `<fasta>.fai`. Other header lines (`@PG`, `@RG`, `@CO`) are passed through.
- **Reference names.** The leading `f`/`r` is stripped from `RNAME` and `RNEXT`, so `fchr1`/`rchr1` both become `chr1`.
- **Read sequence.** `SEQ` (currently the converted read) is replaced with the original from the `YS:Z:` tag:
  - if the read mapped to the reverse strand (flag `0x10`), the original is reverse-complemented to match BWA's orientation;
  - for hard-clipped records (supplementary alignments), the original is trimmed by the leading/trailing `H` lengths in the CIGAR so it matches the stored `SEQ` length;
  - the `YS` tag is then removed and all other tags (`NM`, `AS`, `MD`, …) are kept.

The result is a standard BAM whose coordinates, chromosome names, and read sequences are all in the original reference space, ready for `deamtools bam2bw`, `bam2fragment`, and `qc`.

### Streaming pipeline and threading

`bwa mem` and `samtools sort` run as concurrent subprocesses connected by pipes. A dedicated **feeder thread** converts and writes reads into BWA's stdin while the main thread reads BWA's stdout, performs the restoration, and writes into `samtools sort`'s stdin; an exception in the feeder is propagated and the exit codes of both processes are checked. `--threads` is split between the two stages (`bwa mem` gets ⌈n/2⌉, `samtools sort` the rest). After sorting, the BAM is indexed with `samtools index`.

## Examples

```bash
# Paired-end, 8 threads
deamtools align \
    --fasta hg38.fa \
    --fastq1 sample_R1.fq.gz \
    --fastq2 sample_R2.fq.gz \
    --threads 8 \
    --out_dir results \
    --out_name sample
# -> results/sample.bam (+ results/sample.bam.bai)

# Single-end with a read group
deamtools align \
    --fasta hg38.fa \
    --fastq1 sample.fq.gz \
    --read_group '@RG\tID:s1\tSM:sample1\tLB:lib1\tPL:ILLUMINA' \
    --out_dir results \
    --out_name sample

# Index built in a custom location -> point align at it with --index
deamtools index --fasta hg38.fa --out_dir indexes --out_name hg38
deamtools align \
    --fasta hg38.fa \
    --index indexes/hg38.deamtools.c2t \
    --fastq1 sample_R1.fq.gz --fastq2 sample_R2.fq.gz \
    --out_dir results --out_name sample
```

## Notes

- `--threads` is divided between `bwa mem` and `samtools sort` (roughly half each).
- Paired FASTQs must contain the same number of reads in the same order; a length mismatch raises an error.
- The output BAM is coordinate-sorted and indexed, so it is immediately usable by the rest of the toolkit.
