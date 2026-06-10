# bam2fragment

Convert a coordinate-sorted BAM file to a per-fragment editing-signal table. Each row is a unique fragment defined by its coordinates and the exact set of editing positions it carries — a compact, single-molecule representation suitable for bulk or single-cell analysis.

## Synopsis

```
deamtools bam2fragment --bam FILE --fasta FILE --output FILE [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bam FILE` | Coordinate-sorted, indexed BAM file (`.bai` required). |
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). |
| `--output FILE` | Output table path. If it ends in `.gz`, the file is written gzip-compressed. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--barcode` | *(off)* | Add a barcode column (10x fragments-style ordering). Fragments without the tag get `.`. |
| `--barcode_tag TAG` | `CB` | BAM tag carrying the cell barcode. |
| `--min_mapq INT` | `20` | Minimum read mapping quality. |
| `--min_baseq INT` | `20` | Minimum base quality for a position to count as an editing event. |
| `--threads INT` | `1` | Threads for parallel per-chromosome processing. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## Output format

Tab-delimited, one row per unique fragment signature:

```
# without --barcode
chrom   start   end   count   pos1|pos2|...

# with --barcode (10x ordering)
chrom   start   end   barcode   count   pos1|pos2|...
```

- **`count`** — the number of reads/pairs producing that exact `(coords [, barcode], edits)` signature.
- **edits column** — a `|`-separated list of 0-based reference positions showing a `C→T` (forward read) or `G→A` (reverse read) deamination event. Fragments with no detected edits emit `.`.

Editing is **strand-aware**: forward reads record `C→T`, reverse reads record `G→A`. For properly-paired reads, the two mates are merged into one fragment (`start` = min of the two read starts, `end` = max of the two read ends) and their editing positions are unioned. In an unpaired BAM each read is treated as a single-end fragment.

Reads flagged unmapped, duplicate, QC-fail, secondary, or supplementary are always excluded.

## Examples

```bash
# Bulk fragment table
deamtools bam2fragment \
    --bam sample.bam \
    --fasta hg38.fa \
    --output sample.fragments.tsv

# Single-cell, gzip-compressed, with 10x cell barcodes
deamtools bam2fragment \
    --bam sample.bam \
    --fasta hg38.fa \
    --barcode --barcode_tag CB \
    --output sample.fragments.tsv.gz
```

## Notes

- The fragment table is the natural input for single-molecule and single-cell analyses (per-fragment edit patterns, barcode-level aggregation).
- Output paths ending in `.gz` are written gzip-compressed automatically.
