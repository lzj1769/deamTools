# bam2bw

Convert a coordinate-sorted BAM file to a BigWig track of per-base C→T deamination counts.

## Synopsis

```
deamtools bam2bw --bam FILE --fasta FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bam FILE` | Coordinate-sorted BAM file. Must be accompanied by an index (`.bai`). |
| `--fasta FILE` | Reference FASTA file used during alignment. Must be indexed with `samtools faidx` (`.fai`). |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. The BigWig is written to `<out_dir>/<out_name>.bw`. |

## Optional arguments

### Input / scope

| Argument | Default | Description |
|---|---|---|
| `--chrom_sizes FILE` | *(BAM header)* | Tab-delimited chromosome sizes file (`chrom\tsize`). When omitted, chromosome names and sizes are read from the BAM header. |
| `--regions FILE` | *(whole genome)* | BED file of genomic regions to restrict processing to. Only reads overlapping these intervals are examined, which can substantially reduce runtime for targeted analyses. Overlapping intervals within the same chromosome are merged automatically to prevent double-counting. |

### Signal

| Argument | Default | Description |
|---|---|---|
| `--extend_size INT` | `0` | Symmetrically extend each detected deamination site by INT base pairs in both directions before writing to the BigWig. A value of 50 means each event at position *p* contributes signal to [*p*−50, *p*+50]. Implemented as a box-kernel convolution, so the signal at a position equals the number of events within `extend_size` bases. |

### Quality filters

| Argument | Default | Description |
|---|---|---|
| `--min_mapq INT` | `20` | Minimum read mapping quality (MAPQ). Reads strictly below this value are skipped entirely. |
| `--min_baseq INT` | `20` | Minimum base quality (phred score) at a candidate position. Individual bases below this value are not counted even if the read passes MAPQ filtering. |

Regardless of these thresholds, the following reads are always excluded:

- Unmapped reads (flag `0x4`)
- PCR/optical duplicates (flag `0x400`)
- QC-failed reads (flag `0x200`)
- Secondary alignments (flag `0x100`)
- Supplementary alignments (flag `0x800`)

### Performance

| Argument | Default | Description |
|---|---|---|
| `--threads INT` | `1` | Number of threads for parallel processing. Each thread handles one chromosome independently. |

### Global option (before the subcommand)

| Argument | Default | Description |
|---|---|---|
| `--log_level LEVEL` | `INFO` | Logging verbosity. One of `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## Input file preparation

Before running `bam2bw`, your BAM and FASTA files must be sorted and indexed.

```bash
# Sort and index BAM
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

# Index reference FASTA
samtools faidx hg38.fa

# (Optional) Generate chromosome sizes file from FASTA index
cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
```

## Examples

### Whole-genome run

```bash
deamtools bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --out_dir results \
    --out_name sample
```

Writes `results/sample.bw`. Uses default MAPQ ≥ 20 and base quality ≥ 20. Chromosome sizes are read from the BAM header.

### Restrict to peaks and use stricter quality filters

```bash
deamtools bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --regions peaks.bed \
    --min_mapq 30 \
    --min_baseq 30 \
    --out_dir results \
    --out_name sample_peaks
```

Only reads overlapping intervals in `peaks.bed` are processed, which is much faster than a whole-genome run when peaks cover a small fraction of the genome.

### Extend each deamination event by 50 bp

```bash
deamtools bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --extend_size 50 \
    --out_dir results \
    --out_name sample_extended
```

Each C→T event contributes signal to a 101-bp window centred on the event position. Useful when the raw per-base signal is too sparse for downstream peak calling.

### Parallel chromosome processing with explicit chromosome sizes

```bash
deamtools bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --chrom_sizes hg38.chrom.sizes \
    --threads 8 \
    --out_dir results \
    --out_name sample
```

### Enable debug logging

```bash
deamtools --log_level DEBUG bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --out_dir results \
    --out_name sample
```

Note that `--log_level` is a global flag and must appear **before** the subcommand name.

## Output

The output is a BigWig file in variable-step format. Only positions with non-zero signal are written, keeping file sizes small. The BigWig is directly compatible with genome browsers (IGV, UCSC) and downstream tools (deepTools, MACS3, etc.).

The value at each position is the number of deamination events observed there. If `--extend_size` is used, the value at position *p* is the number of raw events within `extend_size` bases of *p*.

## Choosing parameters

**`--min_mapq`** — A value of 20 (default) retains reads with ≥ 99% mapping accuracy. For repetitive regions or multi-mapping reads, raise to 30. Setting to 0 disables MAPQ filtering.

**`--min_baseq`** — A value of 20 (default) corresponds to 99% base-call accuracy. Lowering increases sensitivity but also increases noise from sequencing errors. Raising above 30 can be overly strict for older sequencing data.

**`--extend_size`** — Set to 0 (default) for the raw single-base signal. For footprinting or broad accessibility analysis, values of 50–200 bp are typical. The optimal value depends on the expected size of accessible regions in your assay.

**`--threads`** — Parallelism is at the chromosome level. Setting `--threads` above the number of chromosomes provides no benefit. For a human genome run, 8–24 threads is a reasonable range.
