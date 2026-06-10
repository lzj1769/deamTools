# qc

Compute quality-control metrics for a deaminase-based chromatin accessibility experiment from a coordinate-sorted BAM and its reference FASTA.

Two files are produced: a machine-readable `<out_dir>/<out_name>.json` and a self-contained, MultiQC-style `<out_dir>/<out_name>.html` report. The HTML embeds the summary figure and documents the meaning of every metric inline.

## Synopsis

```
deamtools qc --bam FILE --fasta FILE --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--bam FILE` | Coordinate-sorted BAM file. Must be accompanied by an index (`.bai`). |
| `--fasta FILE` | Reference FASTA file used during alignment. Must be indexed with `samtools faidx` (`.fai`). |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name (without extension) for the outputs. Writes `<out_dir>/<out_name>.json` and `<out_dir>/<out_name>.html`. |

## Optional arguments

### TSS enrichment

| Argument | Default | Description |
|---|---|---|
| `--tss FILE` | *(disabled)* | BED file of transcription start sites. When supplied, an ATAC-style TSS enrichment score and profile are computed. The TSS is taken as the midpoint of each `(chrom, start, end)` interval. |
| `--tss_flank INT` | `2000` | Half-width in base pairs of the window around each TSS. The profile spans `2 * tss_flank + 1` positions. |

### Quality filters

| Argument | Default | Description |
|---|---|---|
| `--min_mapq INT` | `20` | Minimum read mapping quality (MAPQ). Reads strictly below this value are skipped entirely. |
| `--min_baseq INT` | `20` | Minimum base quality (phred score) at a candidate position. Bases below this value count neither as an editing opportunity nor as an edit. |

Regardless of these thresholds, the following reads are always excluded from the editing and fragment-length metrics:

- Unmapped reads (flag `0x4`)
- PCR/optical duplicates (flag `0x400`)
- QC-failed reads (flag `0x200`)
- Secondary alignments (flag `0x100`)
- Supplementary alignments (flag `0x800`)

The read-count metrics (`total`, `duplicate`, `secondary`, `supplementary`, `unmapped`) report counts *before* filtering, so you can see what fraction of the library was discarded.

### Output control

| Argument | Default | Description |
|---|---|---|
| `--no_plot` | *(off)* | Skip rendering and embedding the summary figure in the HTML report. The JSON and the HTML (tables and descriptions) are still produced. |

### Performance

| Argument | Default | Description |
|---|---|---|
| `--threads INT` | `1` | Number of threads for parallel processing. Each thread handles one chromosome independently. |

### Global option (before the subcommand)

| Argument | Default | Description |
|---|---|---|
| `--log_level LEVEL` | `INFO` | Logging verbosity. One of `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## Input file preparation

```bash
# Sort and index BAM
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

# Index reference FASTA
samtools faidx hg38.fa
```

## Examples

### Core metrics from BAM + FASTA

```bash
deamtools qc \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --out_dir results \
    --out_name sample
```

Produces `results/sample.json` and `results/sample.html`.

### Add TSS enrichment and run on multiple threads

```bash
deamtools qc \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --tss tss.bed \
    --threads 8 \
    --out_dir results \
    --out_name sample
```

### Skip the figure, stricter quality filters

```bash
deamtools qc \
    --bam sample.sorted.bam \
    --fasta hg38.fa \
    --min_mapq 30 \
    --min_baseq 30 \
    --no_plot \
    --out_dir results \
    --out_name sample
```

## Metrics

### Read statistics (`reads`)

Counts of `total`, `passing`, `unmapped`, `duplicate`, `secondary`, `supplementary`, and `proper_pair` reads, plus `duplicate_rate` (over total reads) and `proper_pair_rate` (over passing reads). A low passing fraction or a high duplicate rate points to library-complexity problems.

### Editing statistics (`editing`)

The core signal-quality metrics:

- **`total_opportunities`** — the number of editable reference C/G positions (covered by passing reads, with both flanking bases present, passing `--min_baseq`). Counted strand-agnostically (matching `bam2bw`): every reference **C** *and* every reference **G** the read covers is an opportunity, regardless of read orientation.
- **`total_edits`** — the number of those positions showing a deamination event: a `C→T` mismatch at a reference C or a `G→A` mismatch at a reference G, regardless of read orientation.
- **`global_edit_rate`** — `total_edits / total_opportunities`. The single most important number: a successful deaminase treatment drives this well above the background sequencing-error rate.
- **`mean_edits_per_read`**, **`median_edits_per_read`** — the per-read editing distribution. Deaminase reads typically carry many edits, in contrast to the two Tn5 insertions of a standard ATAC read.

### Per-read edit rate (`edit_rate_per_read`)

The fraction of editable bases that were actually edited, computed **per read**. For each read, the *editable* bases are the reference cytosines and guanines it covers (counted strand-agnostically, gated by `--min_baseq`); the *edited* bases are those showing a `C→T` or `G→A` deamination event. The per-read rate is `edited / editable`.

| Field | Description |
|---|---|
| `n_reads` | Number of reads with at least one editable base (the rest cannot have a rate). |
| `mean`, `median` | Centre of the per-read edit-rate distribution. `mean` is exact; `median` is taken from the histogram bin centres. |
| `histogram` | Counts across 50 equal-width bins spanning the `[0, 1]` rate range. |
| `bin_edges` | The 51 bin boundaries, so `histogram[i]` covers `[bin_edges[i], bin_edges[i+1])`. |

This complements `mean_edits_per_read`: the raw count scales with read length and coverage of editable bases, whereas the rate normalises by how many editable bases each read actually had, making it directly comparable across reads and libraries. A higher, well-separated distribution indicates stronger, more uniform deaminase activity. The distribution is drawn as its own panel in the PNG summary.

### Trinucleotide context bias (`context`)

For each cytosine-centred trinucleotide (e.g. `TCG`, `ACA`), the number of `edits`, `opportunities`, and the resulting `edit_fraction`. `G→A` events are reverse-complemented into the unified `C→T` orientation, so both strands are reported together. This is the enzyme's **sequence-preference fingerprint** — for example, DddA strongly prefers `TC` contexts, while relaxed-bias enzymes such as DddSs/SsdAtox edit more uniformly across contexts. A strongly skewed profile means downstream footprinting will benefit from enzyme-bias correction.

### Fragment-length distribution (`fragment_length`)

`mean`, `median`, and `n_pairs`, computed from `abs(template_length)` of properly-paired read 1 (so each pair is counted once). For an ATAC-style library this should show the characteristic nucleosome-laddering periodicity in the PNG panel.

### TSS enrichment (`tss_enrichment`, optional)

Present only when `--tss` is supplied. The Tn5 insertion 5′ end (`reference_start` on forward reads, `reference_end − 1` on reverse reads) is aggregated into a profile centred on every TSS, normalised by the mean insertion density in the outermost 100 bp of each flank. The `score` is the mean of the normalised profile in the central ±50 bp window; the `profile` array is the full normalised curve. Higher is better — a value above ~6–10 indicates good accessible-chromatin enrichment.

## Output

Two files are written to `--out_dir`:

**`<out_name>.json`** — a machine-readable document with all the sections described above. Suitable for aggregating across many samples (for example, feeding into a comparison table).

**`<out_name>.html`** — a self-contained, MultiQC-style report (no external files or network needed). It opens with headline summary cards, embeds the multi-panel summary figure, and presents every metric in a table alongside a plain-language description of its meaning. The embedded figure (omitted with `--no_plot`) has the panels:

1. Trinucleotide context edit fraction (enzyme fingerprint)
2. Edits-per-read histogram (raw count)
3. Per-read edit-rate distribution (edited / editable)
4. Fragment-length distribution
5. TSS enrichment profile (only when `--tss` is supplied)

## Choosing parameters

**`--min_mapq` / `--min_baseq`** — Keep these consistent with the values used in `bam2bw` / `bam2fragment` so the QC reflects the data your downstream analysis actually sees. Defaults of 20 correspond to ~99% accuracy.

**`--tss_flank`** — 2000 bp (default) matches the conventional ATAC TSS-enrichment window. There is rarely a reason to change it.

**`--threads`** — Parallelism is at the chromosome level; setting `--threads` above the number of chromosomes provides no benefit. The optional TSS-enrichment pass runs separately and is not parallelised.
