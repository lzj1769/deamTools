# footprint

Score transcription-factor footprints at a set of motif sites using a per-base editing BigWig. A bound factor shields its motif from deaminase editing, so a footprint appears as a local depletion of signal at the motif relative to its flanks.

## Synopsis

```
deamtools footprint --bigwig FILE --regions FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bigwig FILE` | Per-base editing BigWig, e.g. produced by [`deamtools bam2bw`](bam2bw.md). |
| `--regions FILE` | BED of motif sites to score (e.g. from [`deamtools match`](match.md)). Column 4, if present, is carried through as the site name. |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. Writes `<out_dir>/<out_name>.bed`. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--n_shuffles INT` | `1000` | Permutations used to build the footprint p-value null (only computed for sites with a positive score). |
| `--threads INT` | `1` | Number of threads; chromosomes are scored in parallel. |
| `--seed INT` | *(unseeded)* | RNG seed for reproducible p-values. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## How it works

For a motif site of width `L = end − start`, the per-base signal is read over the `3 × L` window `[start − L, end + L)` and split into three equal parts:

```
[ left flank ][ motif centre ][ right flank ]
  L bases        L bases          L bases
```

The footprint score is

```
fp_score = mean(left flank) + mean(right flank) − mean(centre)
```

so a depleted centre flanked by high signal yields a positive score. For sites with a positive score, the per-base values in the window are permuted `--n_shuffles` times to form a null distribution, and the p-value is `(#permutations with score ≥ observed + 1) / (n_shuffles + 1)`. Sites with a non-positive score are assigned `p_value = 1.0` without permutation.

Sites whose `3 × L` window would extend beyond the chromosome ends (or are absent from the BigWig header) are skipped.

## Output

A BED-like, tab-delimited file with one line per scored site:

```
chrom    start    end    name    fp_score    p_value
```

- **`fp_score`** — the footprint score (higher = stronger depletion at the motif).
- **`p_value`** — permutation p-value (small = significant footprint).

## Examples

```bash
# Score motif sites from `match` against a bam2bw track
deamtools footprint \
    --bigwig sample.bw \
    --regions mpbs.bed \
    --out_dir results --out_name footprints

# Reproducible p-values, 8 threads
deamtools footprint \
    --bigwig sample.bw \
    --regions mpbs.bed \
    --seed 0 --threads 8 \
    --out_dir results --out_name footprints
```

## Notes

- Use a **count**-mode BigWig (raw or normalized); the score is scale-aware, so normalized tracks are fine for comparing across samples.
- Filtering the output by `p_value` (e.g. `< 0.001`) gives the set of confidently bound sites.
