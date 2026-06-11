# match

Scan the reference sequence of a set of genomic regions (e.g. accessible peaks) for transcription-factor motif occurrences and write them as a BED of motif-predicted binding sites (MPBSs). These sites are the anchors for downstream footprinting and occupancy analysis.

## Synopsis

```
deamtools match --fasta FILE --regions FILE --output FILE [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). |
| `--regions FILE` | BED file of regions to scan. Overlapping intervals are merged. |
| `--output FILE` | Output BED path. Parent directories are created automatically. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--jaspar_release STR` | `JASPAR2024` | JASPAR release to fetch motifs from. |
| `--collection STR` | `CORE` | JASPAR motif collection (e.g. `CORE`, `UNVALIDATED`). |
| `--tax_group GROUP...` | `vertebrates` | One or more JASPAR taxonomic groups. |
| `--p_value FLOAT` | `1e-4` | Significance threshold for motif hits. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## Requirements

Motifs are fetched from JASPAR via the optional **`pyjaspar`** package (`pip install pyjaspar`). MOODS (a core dependency) performs the scanning.

## How it works

Scanning is performed with MOODS:

1. Each motif's count matrix is converted to a **log-odds matrix** against a flat background (with a small pseudocount), and its **reverse complement** is added so both strands are scanned.
2. A per-motif **score threshold** is derived from `--p_value` with `MOODS.tools.threshold_from_p`.
3. For each region, the reference sequence is read from the FASTA and scanned with MOODS; every hit at or above the threshold is reported.

For a hit at sequence position *p* (0-based) of a motif of width *w* in a region starting at genomic coordinate *s*, the reported interval is `[s + p, s + p + w)`, with strand `+` for the forward matrix and `-` for the reverse complement.

## Output

A 6-column BED file, one line per hit:

```
chrom    start    end    motif    score    strand
```

- **`motif`** — the motif label (`<matrix_id>.<name>`, e.g. `MA0139.1.CTCF`).
- **`score`** — the MOODS log-odds bitscore of the match (higher is a better match).
- **`strand`** — `+` (forward matrix) or `-` (reverse complement).

## Examples

```bash
# Scan peaks against JASPAR CORE vertebrate motifs
deamtools match \
    --fasta hg38.fa \
    --regions peaks.bed \
    --output mpbs.bed

# Stricter threshold, explicit collection / taxonomic group
deamtools match \
    --fasta hg38.fa \
    --regions peaks.bed \
    --collection CORE \
    --tax_group vertebrates \
    --p_value 1e-5 \
    --output mpbs.bed
```

## Notes

- Restricting `--regions` to accessible peaks keeps the output to plausible binding sites and is much faster than scanning the whole genome.
- The BED `score` column holds the raw MOODS bitscore (not rescaled to 0–1000), which downstream footprinting can use directly.
