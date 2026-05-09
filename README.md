# DeamTools

A Python command-line toolkit for deaminase-based chromatin accessibility analysis.

DeamTools quantifies cytosine deamination (C→T) events in aligned sequencing reads and converts them into BigWig coverage tracks for genome-browser visualisation and downstream signal analysis. It is designed for deaminase-based assays — such as CHEESE-seq, DamID variants, or any method that uses a deaminase to mark accessible chromatin — where single-base C→T editing events serve as the accessibility signal.

## How it works

For each aligned read, DeamTools scans reference cytosine positions covered by the read and records positions where a C→T (forward strand) or G→A (reverse strand) conversion is observed. These per-base editing counts are accumulated genome-wide and written as a BigWig file.

## Installation

**Requirements:** Python ≥ 3.10, `samtools` (for indexing BAM/FASTA files).

```bash
git clone https://github.com/lzj1769/deamTools.git
cd deamTools
pip install .
```

For development (adds pytest, ruff, black, mypy):

```bash
pip install -e ".[dev]"
```

## Quick start

```bash
# 1. Sort and index your BAM file
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

# 2. Index the reference FASTA
samtools faidx hg38.fa

# 3. Convert to BigWig
deamtools bam2bw \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --output sample.bw
```

## Commands

### Global options

```
deamtools [--version] [--log_level LEVEL] <command>
```

| Option | Default | Description |
|---|---|---|
| `--version` | — | Print version and exit |
| `--log_level` | `INFO` | Verbosity: `DEBUG`, `INFO`, `WARNING`, or `ERROR` |

---

### `bam2bw` — BAM to BigWig

Convert a coordinate-sorted BAM file to a per-base BigWig track of C→T deamination counts.

```
deamtools bam2bw --bam FILE --fasta FILE --output FILE [options]
```

#### Required arguments

| Argument | Description |
|---|---|
| `--bam FILE` | Coordinate-sorted BAM file. Must be indexed (`.bai`). |
| `--fasta FILE` | Reference FASTA file. Must be indexed with `samtools faidx` (`.fai`). |
| `--output FILE` | Output BigWig path (`.bw`). Parent directories are created automatically. |

#### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--chrom_sizes FILE` | *(from BAM header)* | Tab-delimited chromosome sizes file (`chrom\tsize`). If omitted, sizes are inferred from the BAM header. |
| `--regions FILE` | *(whole genome)* | BED file of regions to restrict analysis to. Overlapping intervals are merged automatically. |
| `--extend_size INT` | `0` | Symmetrically extend each detected editing site by INT base pairs before writing to the BigWig. Useful for smoothing sparse signals or defining accessible windows around editing events. |
| `--min_mapq INT` | `20` | Minimum read mapping quality. Reads below this threshold are skipped. |
| `--min_baseq INT` | `20` | Minimum base quality at a position. Bases below this threshold are not counted. |
| `--threads INT` | `1` | Number of threads for parallel chromosome processing. |

Reads that are unmapped, duplicate, QC-failed, secondary, or supplementary are always excluded regardless of quality thresholds.

#### Examples

```bash
# Whole-genome run with default quality thresholds
deamtools bam2bw \
    --bam sample.bam \
    --fasta hg38.fa \
    --output sample.bw

# Restrict to peaks, use stricter filters, run on 4 threads
deamtools bam2bw \
    --bam sample.bam \
    --fasta hg38.fa \
    --regions peaks.bed \
    --min_mapq 30 \
    --min_baseq 30 \
    --threads 4 \
    --output sample_peaks.bw

# Extend each editing site by 50 bp in both directions
deamtools bam2bw \
    --bam sample.bam \
    --fasta hg38.fa \
    --extend_size 50 \
    --output sample_extended.bw

# Provide an explicit chromosome sizes file and enable debug logging
deamtools --log_level DEBUG bam2bw \
    --bam sample.bam \
    --fasta hg38.fa \
    --chrom_sizes hg38.chrom.sizes \
    --output sample.bw
```

## Running tests

```bash
pytest                    # run all tests
pytest -v                 # verbose output
pytest tests/test_bam2bw.py  # run a specific test file
```

## License

MIT — see [LICENSE](LICENSE).
