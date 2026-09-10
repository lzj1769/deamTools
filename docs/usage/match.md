# match

Scan the reference sequence of a set of genomic regions (e.g. accessible peaks) for transcription-factor motif occurrences and write them as a BED of motif-predicted binding sites (MPBSs). These sites are the anchors for downstream footprinting and occupancy analysis.

## Synopsis

```
deamtools match --fasta FILE --regions FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). |
| `--regions FILE` | BED file of regions to scan. Overlapping intervals are merged. |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the output. Writes `<out_dir>/<out_name>.bed`. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--jaspar_release STR` | `JASPAR2024` | JASPAR release to fetch motifs from. |
| `--collection STR` | `CORE` | JASPAR motif collection (e.g. `CORE`, `UNVALIDATED`). |
| `--tax_group GROUP...` | `vertebrates` | One or more JASPAR taxonomic groups. |
| `--p_value FLOAT` | `1e-4` | Significance threshold for motif hits. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## Requirements

Motifs come from one of two places. Pass `--motif_files` to scan with local
motif files, or omit it to fetch from JASPAR via the optional **`pyjaspar`**
package (`pip install pyjaspar`). [motifmatchpy](https://github.com/lzj1769/motifmatchpy)
(a core dependency) performs the scanning.

### Using local motif files

```bash
deamtools match \
    --fasta hg38.fa \
    --regions peaks.bed \
    --out_dir results \
    --out_name mpbs \
    --motif_files motifs/MA0139.1.pfm motifs/MA0095.2.pfm
```

`.pfm` files are position frequency matrices (four whitespace-separated rows,
`A`/`C`/`G`/`T`); `.adm` files are adjacent dinucleotide models. Each motif is
named after its file stem, so `MA0139.1.pfm` appears as `MA0139.1` in the BED
`name` column. `--motif_files` takes precedence over the JASPAR options and
needs no `pyjaspar` install.

## How it works

Scanning is performed with motifmatchpy:

1. Each motif's count matrix is converted to a **log-odds matrix** against a flat background (with a small pseudocount), and its **reverse complement** is added so both strands are scanned.
2. A per-motif **score threshold** is derived from `--p_value` with `motifmatchpy.tools.threshold_from_p`.
3. For each region, the reference sequence is read from the FASTA and scanned on both strands; every hit at or above the threshold is reported.

For a hit at sequence position *p* (0-based) of a motif of width *w* in a region starting at genomic coordinate *s*, the reported interval is `[s + p, s + p + w)`, with strand `+` for the forward matrix and `-` for the reverse complement.

## Output

A 6-column BED file, one line per hit:

```
chrom    start    end    motif    score    strand
```

- **`motif`** — the motif label (`<matrix_id>.<name>`, e.g. `MA0139.1.CTCF`).
- **`score`** — the log-odds bitscore of the match (higher is a better match).
- **`strand`** — `+` (forward matrix) or `-` (reverse complement).

## Examples

```bash
# Scan peaks against JASPAR CORE vertebrate motifs
deamtools match \
    --fasta hg38.fa \
    --regions peaks.bed \
    --out_dir results --out_name mpbs

# Stricter threshold, explicit collection / taxonomic group
deamtools match \
    --fasta hg38.fa \
    --regions peaks.bed \
    --collection CORE \
    --tax_group vertebrates \
    --p_value 1e-5 \
    --out_dir results --out_name mpbs
```

## Notes

- Restricting `--regions` to accessible peaks keeps the output to plausible binding sites and is much faster than scanning the whole genome.
- The BED `score` column holds the raw log-odds bitscore (not rescaled to 0–1000), which downstream footprinting can use directly.
