# index

Build a deamination-aware reference index for a FASTA file. Run this once per reference before `deamtools align`.

`index` produces two things next to the FASTA:

1. the standard FASTA index (`<fasta>.fai`, via `samtools faidx`), used by `bam2bw`, `bam2fragment`, `qc`, and `plot_motif`; and
2. a **doubly-converted BWA index** used by `align` to map heavily deaminated reads.

## Synopsis

```
deamtools index --fasta FILE [--out_dir DIR] [--out_name NAME] [--force]
```

## Arguments


| Argument            | Default                   | Description                                                                                                    |
| ------------------- | ------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `--fasta FILE`      | *(required)*              | Reference FASTA to index.                                                                                      |
| `--out_dir DIR`     | *(the FASTA's directory)* | Directory for the converted reference and BWA index.                                                           |
| `--out_name NAME`   | *(the FASTA file name)*   | Base name for the converted reference and BWA index.                                                           |
| `--force`           | *(off)*                   | Rebuild every output even if it already exists. By default, steps whose output is already present are skipped. |
| `--log_level LEVEL` | `INFO`                    | Global flag (before the subcommand):`DEBUG`, `INFO`, `WARNING`, `ERROR`.                                       |

:::{note}
`--out_dir`/`--out_name` control only the deamtools-specific converted reference and its BWA index. The standard `<fasta>.fai` is **always** written next to the FASTA, because the pysam-based subcommands (`bam2bw`, `bam2fragment`, `qc`, `plot_motif`) require it there. By default the converted index is also written next to the FASTA — which is where `deamtools align` looks for it.
:::

## Requirements

`bwa` and `samtools` must be on your `PATH`.

## How it works

Deaminated reads carry many `C→T` (top strand) or `G→A` (bottom strand) conversions, so they align poorly against an unmodified reference. Following the bwa-meth strategy, `index` reduces the alphabet by writing a *doubly-converted* copy of the reference in which every chromosome appears twice:

- `f<chrom>` — the forward sequence with all **C converted to T**;
- `r<chrom>` — the forward sequence with all **G converted to A**.

`bwa index` is then built on this converted reference. During `align`, read 1 is `C→T`-converted and read 2 is `G→A`-converted, so each read pair maps to a single converted contig (`f…` for top-strand-derived fragments, `r…` for bottom-strand-derived), which preserves proper pairing. The `f`/`r` prefix is stripped from the chromosome name when the final BAM is written.

## Outputs


| File                                                      | Location          | Description                              |
| --------------------------------------------------------- | ----------------- | ---------------------------------------- |
| `<fasta>.fai`                                             | next to the FASTA | Standard FASTA index (`samtools faidx`). |
| `<out_dir>/<out_name>.deamtools.c2t`                      | `--out_dir`       | The doubly-converted reference.          |
| `<out_dir>/<out_name>.deamtools.c2t.{amb,ann,bwt,pac,sa}` | `--out_dir`       | BWA-MEM index files.                     |

With the defaults, `<out_dir>/<out_name>` resolves to the FASTA's own path, so the converted index is written as `<fasta>.deamtools.c2t*` next to the FASTA.

## Examples

```bash
# Build the index next to the FASTA (skips outputs that already exist)
deamtools index --fasta hg38.fa

# Write the converted index to a separate directory
deamtools index --fasta hg38.fa --out_dir indexes --out_name hg38

# Force a full rebuild
deamtools index --fasta hg38.fa --force
```

After this completes, align reads with [`deamtools align`](align.md).

## Notes

- The build is idempotent: re-running without `--force` skips the `faidx`, conversion, and `bwa index` steps whose outputs already exist. Use `--force` after changing the FASTA.
- The converted reference doubles the genome size on disk, and `bwa index` can take a while and use substantial memory for large genomes — this is a one-time cost per reference.
