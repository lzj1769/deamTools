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
```

Reads are converted on the fly to match the doubly-converted reference: read 1 is `C→T`-converted and read 2 (if paired) is `G→A`-converted. The original, unconverted read sequence is stashed in a `YS:Z:` tag and carried through `bwa mem -C`. After alignment, deamtools:

- restores the original SEQ from the `YS` tag (reverse-complemented for reverse-strand reads, trimmed for hard-clipped records) and drops the `YS` tag;
- strips the `f`/`r` prefix from `RNAME`/`RNEXT` so chromosome names match the original reference;
- rewrites the header `@SQ` lines from the original `.fai`.

The result is a standard BAM whose coordinates and sequences are in the original reference space, ready for `deamtools bam2bw`, `bam2fragment`, and `qc`.

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
