# bam2fragment

Convert a coordinate-sorted BAM file to a per-fragment editing-signal table. Each row is a unique fragment defined by its coordinates and the exact set of editing positions it carries — a compact, single-molecule representation suitable for bulk or single-cell analysis.

## Synopsis

```
deamtools bam2fragment --bam FILE --fasta FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bam FILE` | Coordinate-sorted, indexed BAM file (`.bai` required). |
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. Writes `<out_dir>/<out_name>.tsv` (or `.tsv.gz` with `--gzip`). |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--gzip` | *(off)* | Write the table gzip-compressed (`<out_name>.tsv.gz`). |
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
chrom   start   end   count   c2t_pos1|c2t_pos2|...   g2a_pos1|g2a_pos2|...

# with --barcode (10x ordering)
chrom   start   end   barcode   count   c2t_positions   g2a_positions
```

- **`count`** — the number of reads/pairs producing that exact `(coords [, barcode], edits)` signature.
- **`C→T` column** — a `|`-separated list of 0-based reference positions where a reference **C** is read as **T**: the **top** strand was deaminated there.
- **`G→A` column** — the same for a reference **G** read as **A**: the **bottom** strand was deaminated.

Either column is `.` when the fragment shows no event of that kind.

Both patterns are called on **every** read, whatever its orientation, and kept in separate columns so the edited strand is preserved. Read orientation is not what decides the pattern: a BAM stores `SEQ` in reference orientation, and library prep fixes the deaminated U into a real T:A pair that both strands carry, so a read of either orientation reports either pattern. A double-stranded deaminase edits both strands, so one fragment routinely carries both. See [Algorithm → Strand convention](../algorithm.md#strand-convention).

For properly-paired reads, the two mates are merged into one fragment (`start` = min of the two read starts, `end` = max of the two read ends). A position both mates cover is reported once; if they disagree the higher base quality wins, and an equal-quality disagreement is dropped as ambiguous. In an unpaired BAM each read is treated as a single-end fragment.

```{note}
Before 2026-09-10 this command recorded `C→T` only on forward reads and `G→A` only on reverse ones. On real ACCESS-ATAC data that discarded about **36%** of editing events, so tables written by an earlier version are not comparable with current output — and they have one editing column rather than two.
```

Reads flagged unmapped, duplicate, QC-fail, secondary, or supplementary are always excluded.

## Examples

```bash
# Bulk fragment table -> results/sample.tsv
deamtools bam2fragment \
    --bam sample.bam \
    --fasta hg38.fa \
    --out_dir results --out_name sample

# Single-cell, gzip-compressed, with 10x cell barcodes -> results/sample.tsv.gz
deamtools bam2fragment \
    --bam sample.bam \
    --fasta hg38.fa \
    --barcode --barcode_tag CB --gzip \
    --out_dir results --out_name sample
```

## Notes

- The fragment table is the natural input for single-molecule and single-cell analyses (per-fragment edit patterns, barcode-level aggregation).
- Pass `--gzip` to write `<out_name>.tsv.gz` instead of `<out_name>.tsv`.
