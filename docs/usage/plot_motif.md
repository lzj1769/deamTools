# plot_motif

Build a sequence logo of the deaminase's flanking-sequence preference from a per-base BigWig of editing counts. This visualises the enzyme's intrinsic bias (e.g. DddA's strong `TC` preference), which is useful for QC and for deciding whether downstream footprinting needs bias correction.

## Synopsis

```
deamtools plot_motif --bigwig FILE --fasta FILE --output FILE [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bigwig FILE` | Per-base BigWig of editing-event counts, as produced by `deamtools bam2bw --mode count --extend_size 0`. |
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). |
| `--output FILE` | Output plot path. The format is inferred from the extension (`.png`, `.pdf`, `.svg`, …). A `.csv` of the bit-score matrix is written next to it. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--regions FILE` | *(whole BigWig)* | BED file restricting analysis to specific intervals. When omitted, every chromosome in the BigWig header is processed. |
| `--window_size INT` | `10` | Window size in bp around each editing site. Must be ≥ 2. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## How it works

For every base with a non-zero count, the surrounding reference window is read from the FASTA and added to a position-weight matrix, weighted by the count. The **centre base** (the editing site itself) is excluded, so the logo reflects the enzyme's *flanking* preference. When the centre base is a `G`, the window is reverse-complemented before accumulation, so `C→T` and `G→A` events are unified in the `C→T` orientation. Counts are converted to information content (bits) and rendered with `logomaker`; the bit-score matrix is also written to a sibling `.csv`.

## Examples

```bash
# Whole-BigWig deaminase motif logo, default 10-bp window
deamtools plot_motif \
    --bigwig sample.bw \
    --fasta hg38.fa \
    --output motif.png

# Restrict to peaks and use a wider window
deamtools plot_motif \
    --bigwig sample.bw \
    --fasta hg38.fa \
    --regions peaks.bed \
    --window_size 20 \
    --output motif.pdf
```

## Notes

- Use a **count-mode, `--extend_size 0`** BigWig. With a non-zero extension each editing event is broadcast to neighbouring bases, which would inflate non-editing positions and skew the logo.
- A `.csv` of the bit-score matrix is written alongside the plot for downstream analysis.
