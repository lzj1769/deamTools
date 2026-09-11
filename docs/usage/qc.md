# qc

Compute quality-control metrics for a deaminase-based chromatin accessibility experiment from a coordinate-sorted BAM and its reference FASTA.

A machine-readable `<out_dir>/<out_name>.json` and a self-contained, MultiQC-style `<out_dir>/<out_name>.html` report are produced. The HTML embeds the summary figure and documents the meaning of every metric inline. Alongside them, every plotted distribution gets a CSV holding the numbers behind the plot, so it can be re-drawn without rerunning: `<out_name>.edits_per_fragment.csv` and `<out_name>.edit_rate_per_fragment.csv` always, and `<out_name>.tss_enrichment.csv` with `--tss`. These are written even with `--no_plot`.

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
| `--out_name NAME` | Base name (without extension) for the outputs. Writes `<out_dir>/<out_name>.json` and `<out_dir>/<out_name>.html`, the two distribution CSVs, and `<out_name>.tss_enrichment.csv` when `--tss` is given. |

## Optional arguments

### TSS enrichment

| Argument | Default | Description |
|---|---|---|
| `--tss FILE` | *(disabled)* | BED file of transcription start sites. When supplied, a TSS enrichment score and aggregate profile are computed. The TSS is the midpoint of each `(chrom, start, end)` interval. Strand is read from column 6, or from column 4 for four-column `chrom start end strand` files; minus-strand TSS are flipped so upstream is always on the left. |
| `--tss_flank INT` | `2000` | Half-width in base pairs of the window around each TSS, rounded down to a whole number of 10 bp bins. The profile spans `2 * tss_flank / 10` bins. |

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

### Subsampling a large BAM

`--n_reads` computes the QC metrics from a random subset instead of every read,
which is worth doing when a full pass is slow:

```bash
deamtools qc --bam sample.bam --fasta hg38.fa \
    --out_dir results --out_name sample --n_reads 5000000
```

The sampling fraction is `n_reads / (reads in the BAM)`, read straight from the
BAM index, so nothing is scanned twice. Each **fragment** is kept independently
with that probability, which makes the sample uniform across the genome rather
than biased toward the first chromosomes. The draw is keyed on the read name, so
both mates of a pair always share their fragment's fate — drawing the two mates
separately would leave most kept fragments with only one mate and roughly halve
every per-fragment edit count. The draw is seeded, so rerunning the same command
samples the same fragments.

**The draw is blind to what a read contains.** The coin is flipped before the
read is examined, so a read carrying no editing event is kept at exactly the
same rate as one full of them — nothing is filtered on edits, mapping quality,
or anything else at sampling time. (The usual filters still apply afterwards,
to the sampled reads, exactly as in a full run.)

Rates and distributions are unbiased under this sampling: the editing rate,
duplicate rate, trinucleotide context bias, motif PWM and fragment-length
distribution all mean the same thing as in a full run. **The absolute counts are
counts of the sample**, not estimates of the whole file; the `sampling` block of
the JSON records `fraction` and `reads_in_bam` so they can be scaled if needed.

TSS enrichment always uses every read in its windows. It is a ratio computed
over a small part of the genome, so subsampling it would add noise without
saving meaningful time.

## Output control

| Argument | Default | Description |
|---|---|---|
| `--no_plot` | *(off)* | Skip rendering and embedding the summary figure in the HTML report. The JSON, the CSVs and the HTML (tables and descriptions) are still produced. |
| `--logo_scale {bits,frequency}` | `bits` | Y axis of the deaminase motif logo. `bits` plots information content with the edited base left out — it is always C and would take the full 2 bits, flattening flanks that rarely reach 0.15. `frequency` plots each base's frequency per offset on a 0–1 axis, with the target C drawn at position 0; it is usually the easier of the two to read for a weakly specific enzyme. |

### Performance

| Argument | Default | Description |
|---|---|---|
| `--threads INT` | `1` | Number of worker processes. Each chromosome is processed independently, so the work spreads across cores. On the 3.2 M-record single-end concurrent ACCESS-ATAC BAM, `--threads 12` takes 33 s against 192 s before the switch from threads to processes. |

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

### Add TSS enrichment and run on several worker processes

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

### Library layout (`library_layout`)

Before the main pass, `qc` reads the first 10,000 records of the BAM and calls the library **`paired-end`** if any of them carries the paired flag (`0x1`), otherwise **`single-end`**. One record is enough either way: a single-end library never sets the flag, and a paired-end one sets it on essentially every record, unmapped reads and orphans included. The layout is shown in the report header and recorded at the top of the JSON.

For a **single-end** library the pair-dependent metrics are **left out** rather than reported as zero — `proper_pair`, `proper_pair_rate`, and the whole `fragment_length` block, along with its panel in the summary figure. A proper-pair rate of 0 would read as a mapping failure, and a fragment length of 0 bp is not a length; neither applies when there are no pairs. Everything else — editing, context, motif, TSS enrichment — is computed identically for both layouts, since each single-end read is simply a fragment of one record.

### Read statistics (`reads`)

Counts of `total`, `passing`, `unmapped`, `duplicate`, `secondary`, and `supplementary` reads plus `duplicate_rate` (over total reads); for a paired-end library also `proper_pair` and `proper_pair_rate` (over passing reads). A low passing fraction or a high duplicate rate points to library-complexity problems.

### Fragments (`fragments`)

**Editing is counted per fragment, not per read.** For paired-end data the two mates
are merged before anything is counted, so a reference position that both mates
cover is one observation instead of two. Where the mates disagree at an overlap,
the higher base quality wins — library prep turns a deaminated C into a real T:A
pair that *both* mates then report, so a disagreement there means one of them
misread. Single-end reads are fragments of one record, so nothing changes for them.

| Field | Description |
|---|---|
| `total` | Fragments the editing metrics were computed over. |
| `from_mate_pairs` | Fragments built by merging two mates. |
| `from_single_records` | Fragments that were a single record — unpaired reads, and reads whose mate was unmapped, filtered out, or on another contig. |

The `reads` block above still counts *records*, so for a paired library
`reads.passing` is roughly twice `fragments.total`.

This matters more than it sounds. On `data/ACCESS-ATAC/chr10.bam`, **24.5%** of
`total_opportunities` were mate overlap being counted twice (62.8 M → 47.4 M), and
removing it moved `global_edit_rate` from 0.0824 to 0.0856 — the overlap is not a
random subset, so double-weighting it biased the rate, it did not merely inflate
the counts.

### Editing statistics (`editing`)

The core signal-quality metrics:

- **`total_opportunities`** — the number of editable reference C/G positions (covered by passing fragments, with both flanking bases present, passing `--min_baseq`). Counted strand-agnostically (matching `bam2bw`): every reference **C** *and* every reference **G** the fragment covers is an opportunity, regardless of read orientation. A position covered by both mates counts once.
- **`total_edits`** — the number of those positions showing a deamination event: a `C→T` mismatch at a reference C or a `G→A` mismatch at a reference G, regardless of read orientation.
- **`global_edit_rate`** — `total_edits / total_opportunities`. The single most important number: a successful deaminase treatment drives this well above the background sequencing-error rate.
- **`mean_edits_per_fragment`**, **`median_edits_per_fragment`** — the per-fragment editing distribution. Deaminase fragments typically carry many edits, in contrast to the two Tn5 insertions of a standard ATAC read.
- **`edits_per_fragment_csv`** — the file holding the full distribution (see [Output](#output)).

### Per-fragment edit rate (`edit_rate_per_fragment`)

The fraction of editable bases that were actually edited, computed **per fragment**. For each fragment, the *editable* bases are the distinct reference cytosines and guanines it covers (counted strand-agnostically, gated by `--min_baseq`); the *edited* bases are those showing a `C→T` or `G→A` deamination event. The rate is `edited / editable`.

| Field | Description |
|---|---|
| `n_fragments_with_editable_bases` | Fragments covering at least one reference C or G — the rest have no denominator and so no rate. This counts editable *bases*, not editing *events*: a fragment with no edit at all still contributes, at rate 0. Unrelated to the `--n_reads` subsampling flag. |
| `mean`, `median` | Centre of the per-fragment edit-rate distribution. `mean` is exact; `median` is taken from the histogram bin centres. |
| `histogram_csv` | The file holding the histogram: 200 equal-width bins over `[0, 1]` (see [Output](#output)). Like the TSS profile, the histogram itself is not repeated in the JSON. |

This complements `mean_edits_per_fragment`: the raw count scales with fragment length and coverage of editable bases, whereas the rate normalises by how many editable bases each fragment actually had, making it directly comparable across fragments and libraries. A higher, well-separated distribution indicates stronger, more uniform deaminase activity. The distribution is drawn as its own panel in the PNG summary, on a log-like x axis (linear below 0.01) because real rates pile up well below 0.1.

### Trinucleotide context bias (`context`)

For each cytosine-centred trinucleotide (e.g. `TCG`, `ACA`), the number of `edits`, `opportunities`, and the resulting `edit_fraction`. `G→A` events are reverse-complemented into the unified `C→T` orientation, so both strands are reported together. This is the enzyme's **sequence-preference fingerprint** — for example, DddA strongly prefers `TC` contexts, while relaxed-bias enzymes such as DddSs/SsdAtox edit more uniformly across contexts. A strongly skewed profile means downstream footprinting will benefit from enzyme-bias correction.

### Fragment-length distribution (`fragment_length`, paired-end only)

`mean`, `median`, and `n_pairs`, computed from `abs(template_length)` of properly-paired read 1 (so each pair is counted once). For an ATAC-style library this should show the characteristic nucleosome-laddering periodicity in the PNG panel.

### TSS enrichment (`tss_enrichment`, optional)

Present only when `--tss` is supplied. Computed the way the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) defines it:

1. Take a ±`--tss_flank` window (ENCODE: 2 kb) around each TSS and bin it at 10 bp.
2. Count insertion 5′ ends — `reference_start` on forward reads, `reference_end − 1` on reverse reads — into those bins, flipping the window for minus-strand TSS so upstream is always on the left.
3. Average over TSS, then apply the Greenleaf normalisation: divide by the mean of the two edge means, each over the outermost 100 bp, so the flanks sit at 1.
4. The `score` is the **peak** of that normalised profile.

| Field | Description |
|---|---|
| `score` | Peak of the normalised profile. Higher is better; ENCODE calls ≥5 acceptable and ≥7 ideal for human ATAC, and the ACCESS-ATAC preprint reports ~13–14 for a good library. |
| `n_tss` | TSS that contributed — those on a contig present in the BAM header whose full window fits inside it. |
| `flank`, `bin_size` | Window half-width and bin width actually used, in bp. |
| `background` | Mean insertions per TSS per bin in the outermost 100 bp on each side; the divisor in step 3. |
| `total_insertions` | Insertion sites counted across all TSS windows. |
| `profile_csv` | Name of the CSV holding the profile itself. |

The profile is not duplicated in the JSON. It lives in **`<out_name>.tss_enrichment.csv`**, one row per bin with columns `position` (bin centre, bp from the TSS), `insertions` (raw count summed over TSS), `mean_insertions_per_tss`, and `normalized` (the plotted curve).

One deviation from ENCODE is deliberate. ENCODE reaches the insertion site indirectly, asking `metaseq` for read *coverage* shifted by `−read_len/2` so each read's interval is centred on its cut site; that spreads every insertion over a read-length-wide box, smoothing the profile and depressing the peak. `deamtools` counts the cut site itself at 1 bp — which is what that shift is approximating — so scores run slightly above ENCODE's for the same library, by more the longer the reads.

## Output

**`<out_name>.json`** — a machine-readable document with all the sections described above. Suitable for aggregating across many samples (for example, feeding into a comparison table).

**`<out_name>.html`** — a self-contained, MultiQC-style report (no external files or network needed). It opens with headline summary cards, embeds the plots, and presents every metric in a table alongside a plain-language description of its meaning. The embedded figures (omitted with `--no_plot`) are the four-panel summary figure:

1. Trinucleotide context edit fraction (enzyme fingerprint)
2. Edits-per-fragment histogram (raw count)
3. Per-fragment edit-rate distribution (edited / editable)
4. Fragment-length distribution (paired-end libraries only)

plus the deaminase sequence-motif logo (in bits or frequency, per `--logo_scale`), and — when `--tss` is supplied — the TSS enrichment profile, with the score marked on the curve.

The numbers behind each plot are written as CSV, whether or not the figures are rendered:

| File | Columns | Notes |
|---|---|---|
| `<out_name>.edits_per_fragment.csv` | `edits`, `fragments`, `fraction`, `is_overflow` | One row per edit count from 0 to 100. The last row is an **overflow bin** counting every fragment with at least 100 edits, flagged `is_overflow = True`. |
| `<out_name>.edit_rate_per_fragment.csv` | `bin_start`, `bin_end`, `fragments`, `fraction` | 200 bins of width 0.005. Half-open `[bin_start, bin_end)`, except the last, which also holds a rate of exactly 1. |
| `<out_name>.motif_pfm.csv` | `position`, `A`, `C`, `G`, `T` | Base **counts** at each offset from the edited base, −5 to +5, in the C→T orientation (G→A events reverse-complemented). Position 0 is the edited base itself, filled in as `C = n_events` because it is always C in that orientation. Either logo can be redrawn from it: `pd.read_csv(path, index_col='position')` is the position × base matrix, and dividing each row by its sum gives the frequency logo. |
| `<out_name>.tss_enrichment.csv` | `position`, `insertions`, `mean_insertions_per_tss`, `normalized` | Only with `--tss`; columns as described under TSS enrichment above. |

The JSON names each file (`editing.edits_per_fragment_csv`, `edit_rate_per_fragment.histogram_csv`, `motif.pfm_csv`, `tss_enrichment.profile_csv`) rather than repeating the numbers, so there is one copy of each distribution.

The report opens with a table of the run's provenance — sample, library layout, the BAM, FASTA and (if given) TSS BED paths, the deamtools version and the time it was generated.

## Choosing parameters

**`--min_mapq` / `--min_baseq`** — Keep these consistent with the values used in `bam2bw` / `bam2fragment` so the QC reflects the data your downstream analysis actually sees. Defaults of 20 correspond to ~99% accuracy.

**`--tss_flank`** — 2000 bp (default) is the ENCODE window. Changing it changes the background, since that is defined relative to the window edges, so a score computed with a different flank is not comparable to a published one.

**`--threads`** — The number of worker processes. Parallelism is at the chromosome level, so setting `--threads` above the number of chromosomes provides no benefit, and the largest chromosome sets the floor on run time. Each worker holds one chromosome's reference sequence, so memory grows with the worker count. The optional TSS-enrichment pass runs separately and is not parallelised.
