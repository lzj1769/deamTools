# align

Align deaminated sequencing reads (single- or paired-end) to a reference prepared with [`deamtools index`](index.md), writing a coordinate-sorted, indexed BAM.

## Synopsis

```
deamtools align --fasta FILE --read1 FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--fasta FILE` | Reference FASTA. Must already have been indexed with `deamtools index` (i.e. `<fasta>.fai` and `<fasta>.deamtools.c2t.*` exist). |
| `--read1 FILE` | FASTQ for read 1, or the only FASTQ for single-end input. Plain or gzipped. |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. Writes a sorted, indexed `<out_dir>/<out_name>.bam` (and `.bam.bai`). |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--read2 FILE` | *(single-end)* | FASTQ for read 2. Provide it for paired-end alignment; omit for single-end. |
| `--index FILE` | *(next to the FASTA)* | Path to the converted reference built by `deamtools index` (`<out_dir>/<out_name>.deamtools.c2t`). Use this when the index was built with a custom `--out_dir`/`--out_name`. |
| `--read_group STR` | *(none)* | Read-group line passed to `bwa mem -R`, e.g. `'@RG\tID:s1\tSM:sample1\tLB:lib1\tPL:ILLUMINA'`. |
| `--threads INT` | `1` | Threads used by `bwa mem` and (separately) `samtools sort`. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

:::{note}
By default `align` looks for the converted index next to the FASTA (`<fasta>.deamtools.c2t`). If you built it with a custom location — `deamtools index --fasta ref.fa --out_dir idx --out_name foo` — point `align` at it with `--index idx/foo.deamtools.c2t`. The `index` command logs this exact path (`Converted index: …`) so you can copy it.
:::

## Requirements

`bwa` and `samtools` must be on your `PATH`, and the reference must have been prepared with `deamtools index` first.

## How it works

```
FASTQ(s) ──[convert ×2]──▶ bwa mem -C ──[group + pick best]──▶ <out_name>.sam
   2 candidates / read       |  three-letter   |  keep higher-scoring candidate
   + stash original (YS)     |  alignment to   |  restore SEQ, strip f/r prefix
   + mark candidate (YC)     |  f/r conv. ref  |  rebuild header
                                                        │
                                       samtools sort ───┴──▶ <out_name>.bam (+ .bai)
```

Every read is emitted in **two** converted forms; bwa maps both, and the better-scoring one is kept. The chosen alignments are restored and written to `<out_dir>/<out_name>.sam`, which `samtools sort` then converts to a coordinate-sorted, indexed BAM. (The intermediate `.sam` is left in place.)

### Why two conversions per read

The ACCESS-ATAC deaminase edits cytosines on **both** strands of accessible DNA. After PCR a single read can therefore show **both** `C→T` (from same-strand C deamination) **and** `G→A` (from the complementary strand's C deamination, read as a G→A relative to the forward reference). Converting a read in only one direction (the bwa-meth assumption) would leave the other direction as dense mismatches, hurting mappability for the most-edited (most accessible) reads.

DeamTools instead aligns each read in **both** conversion directions and keeps whichever maps better, so a read dominated by either edit type is recovered. The reference (built by [`deamtools index`](index.md)) is doubly converted: `f<chrom>` = `C→T` of the forward sequence, `r<chrom>` = `G→A` of the forward sequence.

### 1. On-the-fly dual conversion

As reads are streamed to the aligner, each is written twice. The two converted copies share the **original read name**, carry the original sequence in a `YS:Z:` tag, and are marked with a `YC:Z:` candidate tag (both passed through `bwa mem -C`):

- **Single-end** — two candidates per read: `C→T` (`YC:Z:ct`) and `G→A` (`YC:Z:ga`).
- **Paired-end** — two *fragment orientations*, with a consistent direction for the whole pair, written interleaved so `bwa mem -p` pairs each:
  - `f` = (read 1 `C→T`, read 2 `G→A`)
  - `r` = (read 1 `G→A`, read 2 `C→T`)

```
FASTQ record written to bwa:                Resulting (pre-restore) SAM tags:
@read1  YS:Z:ACGTTCGA  YC:Z:ct              SEQ = ATGTTTGA   (C→T converted)
ATGTTTGA                                    YS:Z:ACGTTCGA    (original, via -C)
+                                           YC:Z:ct          (candidate marker)
IIIIIIII
```

Reads with no quality string (e.g. FASTA input) get a placeholder quality.

### 2. Mapping and picking the best candidate

`bwa mem -C -t <n> [-p] [-R <rg>] <converted_ref> -` maps all candidates against the doubly-converted reference. bwa-mem preserves input order, so a read's two candidates (which share the read name) come out consecutively. Records are grouped by read name and partitioned by their `YC` tag; the candidate with the higher **primary alignment score** is kept — for pairs that is the sum of the two mates' `AS`, so a single consistent orientation (`f` or `r`) is chosen for the whole fragment. The losing candidate's records are dropped. (On a tie the first orientation, `f`/`ct`, wins.)

### 3. Restoring the alignments

The kept records are rewritten back into the original reference space:

- **Header.** BWA's `@SQ`/`@HD` lines describe the doubled `f`/`r` contigs, so they are dropped and replaced with a fresh `@HD VN:1.6 SO:coordinate` plus one `@SQ` per chromosome read from the original `<fasta>.fai`. Other header lines (`@PG`, `@RG`, `@CO`) are passed through.
- **Reference names.** The leading `f`/`r` is stripped from `RNAME` and `RNEXT`, so `fchr1`/`rchr1` both become `chr1`.
- **Read sequence.** `SEQ` (the converted read) is replaced with the original from the `YS:Z:` tag — reverse-complemented for reverse-strand records (flag `0x10`), and trimmed by the CIGAR hard-clip lengths for hard-clipped records. The `YS` and `YC` tags are removed; all other tags (`NM`, `AS`, `MD`, …) are kept.

The result is a standard BAM whose coordinates, chromosome names, and read sequences are all in the original reference space, ready for `deamtools bam2bw`, `bam2fragment`, and `qc`.

### Pipeline and threading

`bwa mem` runs as a subprocess: a dedicated **feeder thread** writes both converted candidates of every read into BWA's stdin while the main thread reads BWA's stdout, groups records by read name to pick the best candidate, restores them, and writes the result to `<out_name>.sam`; an exception in the feeder is propagated and BWA's exit code is checked. The SAM is then converted to a coordinate-sorted BAM with `samtools sort` and indexed with `samtools index`. Both `bwa mem` and `samtools sort` use the full `--threads` count (they run one after the other). Emitting two candidates per read roughly doubles the bwa input.

## Examples

```bash
# Paired-end, 8 threads
deamtools align \
    --fasta hg38.fa \
    --read1 sample_R1.fq.gz \
    --read2 sample_R2.fq.gz \
    --threads 8 \
    --out_dir results \
    --out_name sample
# -> results/sample.bam (+ results/sample.bam.bai)

# Single-end with a read group
deamtools align \
    --fasta hg38.fa \
    --read1 sample.fq.gz \
    --read_group '@RG\tID:s1\tSM:sample1\tLB:lib1\tPL:ILLUMINA' \
    --out_dir results \
    --out_name sample

# Index built in a custom location -> point align at it with --index
deamtools index --fasta hg38.fa --out_dir indexes --out_name hg38
deamtools align \
    --fasta hg38.fa \
    --index indexes/hg38.deamtools.c2t \
    --read1 sample_R1.fq.gz --read2 sample_R2.fq.gz \
    --out_dir results --out_name sample
```

## Notes

- `bwa mem` and `samtools sort` each use the full `--threads` count (they run one after the other, not concurrently).
- The intermediate `<out_name>.sam` is kept alongside the BAM; delete it once you have the sorted BAM if you don't need it.
- Paired FASTQs must contain the same number of reads in the same order; a length mismatch raises an error.
- The output BAM is coordinate-sorted and indexed, so it is immediately usable by the rest of the toolkit.
