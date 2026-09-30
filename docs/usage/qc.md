# qc

Compute quality-control metrics for a deaminase-based chromatin accessibility experiment from a coordinate-sorted BAM and its reference FASTA.

A machine-readable `<out_dir>/<out_name>.json` and a self-contained, MultiQC-style `<out_dir>/<out_name>.html` report are produced. The HTML is an offline dashboard with an overall data-check summary, compact KPI cards, sticky navigation, independent interactive SVG panels, and collapsible metric definitions. Alongside them, every plotted distribution gets a CSV holding the numbers behind the plot, so it can be re-drawn without rerunning: `<out_name>.edits_per_fragment.csv` and `<out_name>.edit_rate_per_fragment.csv` always, and `<out_name>.tss_enrichment.csv` with `--tss`. These are written even with `--no_plot`.

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
| `--no_plot` | *(off)* | Skip all plots in the HTML report. The JSON, the CSVs and the HTML (tables and descriptions) are still produced. |
| `--logo_scale {frequency,bits}` | `frequency` | Y axis of the deaminase motif logo. `frequency` plots each base's frequency per offset on a 0–1 axis, with the target C drawn at position 0. `bits` plots information content with the edited base left out — it is always C and would take the full 2 bits. Bits is the more conventional logo, but a deaminase's flanking preference is weak (the flanks rarely reach 0.15 bits), so its letters come out barely legible; that is why `frequency` is the default. |

| `--report-dir DIR` | *(disabled)* | Write `DIR/report.html`, external PNG/SVG assets, companion CSVs and `metrics.json`. JSON/CSV outputs in `--out_dir` remain available; no separate standalone HTML is written there in this mode. |

| `--redact-paths` | *(off)* | Remove directory components from report provenance and recorded command arguments. Applies to HTML and portable `metrics.json`; local QC JSON retains full provenance. |

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

For a **single-end** library the pair-dependent metrics are **left out** rather than reported as zero — `proper_pair`, `proper_pair_rate`, and the whole `fragment_length` block, along with its plot. A proper-pair rate of 0 would read as a mapping failure, and a fragment length of 0 bp is not a length; neither applies when there are no pairs. Everything else — editing, context, motif, TSS enrichment — is computed identically for both layouts, since each single-end read is simply a fragment of one record.

### Read statistics (`reads`)

Sampled counts of `total`, `passing`, `unmapped`, `duplicate`, `secondary`, `supplementary`, `qcfail`, and `low_mapq` reads plus `duplicate_rate` (over total reads); for a paired-end library also `proper_pair` and `proper_pair_rate` (over passing reads). A low passing fraction or a high duplicate rate points to library-complexity problems.

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

- **`total_opportunities`** — the number of editable reference C/G positions (covered by passing fragments, with an A/C/G/T read base passing `--min_baseq`, including contig boundaries). Counted strand-agnostically (matching `bam2bw`): every reference **C** *and* every reference **G** the fragment covers is an opportunity, regardless of read orientation. A position covered by both mates counts once.
- **`total_edits`** — the number of those positions showing a deamination event: a `C→T` mismatch at a reference C or a `G→A` mismatch at a reference G, regardless of read orientation.
- **`global_edit_rate`** — `total_edits / total_opportunities`. Overall fraction of covered editable C/G positions showing C→T or G→A editing. Interpret relative to matched controls and library conditions.
- **`mean_edits_per_fragment`**, **`median_edits_per_fragment`** — the per-fragment editing distribution. Deaminase fragments typically carry many edits, in contrast to the two Tn5 insertions of a standard ATAC read.
- **`edits_per_fragment_csv`** — the file holding the full distribution (see [Output](#output)).

### Per-fragment edit rate (`edit_rate_per_fragment`)

The fraction of editable bases that were actually edited, computed **per fragment**. For each fragment, the *editable* bases are the distinct reference cytosines and guanines it covers (counted strand-agnostically, gated by `--min_baseq`); the *edited* bases are those showing a `C→T` or `G→A` deamination event. The rate is `edited / editable`.

| Field | Description |
|---|---|
| `n_fragments_with_editable_bases` | Fragments covering at least one reference C or G — the rest have no denominator and so no rate. This counts editable *bases*, not editing *events*: a fragment with no edit at all still contributes, at rate 0. Unrelated to the `--n_reads` subsampling flag. |
| `mean`, `median` | Centre of the per-fragment edit-rate distribution. `mean` is exact; `median` is taken from the histogram bin centres. |
| `histogram_csv` | The file holding the histogram: 200 equal-width bins over `[0, 1]` (see [Output](#output)). Like the TSS profile, the histogram itself is not repeated in the JSON. |

This complements `mean_edits_per_fragment`: the raw count scales with fragment length and coverage of editable bases, whereas the rate normalises by how many editable bases each fragment actually had, making it directly comparable across fragments and libraries. A higher, well-separated distribution indicates stronger, more uniform deaminase activity. The distribution is drawn in an independent SVG panel, on a log-like x axis (linear below 0.01) because real rates pile up well below 0.1.

### Trinucleotide context bias (`context`)

For each cytosine-centred trinucleotide (e.g. `TCG`, `ACA`), the number of `edits`, `opportunities`, and the resulting `edit_fraction`. `G→A` events are reverse-complemented into the unified `C→T` orientation, so both strands are reported together. This is the enzyme's **sequence-preference fingerprint** — for example, DddA strongly prefers `TC` contexts, while relaxed-bias enzymes such as DddSs/SsdAtox edit more uniformly across contexts. A strongly skewed profile means downstream footprinting will benefit from enzyme-bias correction.

### Fragment-length distribution (`fragment_length`, paired-end only)

`mean`, `median`, and `n_pairs`, computed from `abs(template_length)` of properly-paired read 1 (so each pair is counted once). For an ATAC-style library this should show the characteristic nucleosome-laddering periodicity in its SVG panel.

### TSS enrichment (`tss_enrichment`, optional)

Present only when `--tss` is supplied. Computed the way the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) defines it:

1. Take a ±`--tss_flank` window (ENCODE: 2 kb) around each TSS and bin it at 10 bp.
2. Count insertion 5′ ends — `reference_start` on forward reads, `reference_end − 1` on reverse reads — into those bins, flipping the window for minus-strand TSS so upstream is always on the left.
3. Average over TSS, then apply the Greenleaf normalisation: divide by the mean of the two edge means, each over the outermost 100 bp, so the flanks sit at 1.
4. The `score` is the **peak** of that normalised profile.

| Field | Description |
|---|---|
| `score` | Peak of the normalised profile. Compare matched protocols and parameters; this implementation does not apply universal pass/fail thresholds. |
| `n_tss` | TSS that contributed — those on a contig present in the BAM header whose full window fits inside it. |
| `flank`, `bin_size` | Window half-width and bin width actually used, in bp. |
| `background` | Mean insertions per TSS per bin in the outermost 100 bp on each side; the divisor in step 3. |
| `total_insertions` | Insertion sites counted across all TSS windows. |
| `profile_csv` | Name of the CSV holding the profile itself. |

The profile is not duplicated in the JSON. It lives in **`<out_name>.tss_enrichment.csv`**, one row per bin with columns `position` (bin centre, bp from the TSS), `insertions` (raw count summed over TSS), `mean_insertions_per_tss`, and `normalized` (the plotted curve).

One deviation from ENCODE is deliberate. ENCODE reaches the insertion site indirectly, asking `metaseq` for read *coverage* shifted by `−read_len/2` so each read's interval is centred on its cut site; that spreads every insertion over a read-length-wide box, smoothing the profile and depressing the peak. `deamtools` counts the cut site itself at 1 bp — which is what that shift is approximating — so scores run slightly above ENCODE's for the same library, by more the longer the reads.

## Output

**`<out_name>.json`** — a machine-readable document with all the sections described above. Suitable for aggregating across many samples (for example, feeding into a comparison table).

**`<out_name>.html`** — an offline dashboard with no CDN or network dependency.
Its default structure is:

1. Compact metadata, Overall QC and six KPI cards (percentages, K/M/B counts,
   precise values and definitions on hover/focus).
2. Editing signal: editable C/G → edits → global rate, with per-fragment summaries
   and separate count/rate distribution panels.
3. Enzyme sequence bias: sequence logo, context rates and a context summary.
4. TSS enrichment: profile, score, background reference and detailed methods.
5. Read retention and alignment diagnostics; fragment length only for paired-end.
6. Run information, complete paths, parameters, recorded command and CSV downloads.
7. Data-driven downstream recommendations.

Overview, editing signal, TSS and the enzyme sequence-bias summary are visible
by default. Detailed motif/context plots and tables remain collapsed. Navigation opens the necessary parent
panels. Each numeric SVG panel supports exact-value hover, horizontal zoom/pan,
reset, SVG export and CSV download. The edit-rate axis is linear through 1%, then
logarithmic. The logo remains a separately downloadable PNG. These features use
small embedded scripts, not a Plotly/CDN dependency. Without JavaScript, charts,
native SVG hover, CSV downloads and expandable tables remain available. Printing
expands the details and hides controls. Layout adapts to mobile screens.

The numbers behind each plot are written as CSV, whether or not the figures are rendered:

| File | Columns | Notes |
|---|---|---|
| `<out_name>.edits_per_fragment.csv` | `edits`, `fragments`, `fraction`, `is_overflow` | One row per edit count from 0 to 100. The last row is an **overflow bin** counting every fragment with at least 100 edits, flagged `is_overflow = True`. |
| `<out_name>.edit_rate_per_fragment.csv` | `bin_start`, `bin_end`, `fragments`, `fraction` | 200 bins of width 0.005. Half-open `[bin_start, bin_end)`, except the last, which also holds a rate of exactly 1. |
| `<out_name>.motif_pfm.csv` | `position`, `A`, `C`, `G`, `T` | Base **counts** at each offset from the edited base, −5 to +5, in the C→T orientation (G→A events reverse-complemented). Position 0 is the edited base itself, filled in as `C = n_events` because it is always C in that orientation. Either logo can be redrawn from it: `pd.read_csv(path, index_col='position')` is the position × base matrix, and dividing each row by its sum gives the frequency logo. |
| `<out_name>.tss_enrichment.csv` | `position`, `insertions`, `mean_insertions_per_tss`, `normalized` | Only with `--tss`; columns as described under TSS enrichment above. |

The JSON names each file (`editing.edits_per_fragment_csv`, `edit_rate_per_fragment.histogram_csv`, `motif.pfm_csv`, `tss_enrichment.profile_csv`) rather than repeating the numbers, so there is one copy of each distribution.

The header shows filenames, sample, layout, recorded version and generation time.
Complete paths are under **Show full paths**. Collapsing is not redaction: paths
remain in the default HTML. Use `--redact-paths` to remove directory components. The CLI records its actual argument vector in provenance;
its shell-quoted representation is displayed as Command. Original shell quoting
cannot be recovered. Library calls omit Command unless it is explicitly supplied.
Missing metadata in older reports is shown as “Not recorded”, never inferred.

## Choosing parameters

**`--min_mapq` / `--min_baseq`** — Keep these consistent with the values used in `bam2bw` / `bam2fragment` so the QC reflects the data your downstream analysis actually sees. Defaults of 20 correspond to ~99% accuracy.

**`--tss_flank`** — 2000 bp (default) is the ENCODE window. Changing it changes the background, since that is defined relative to the window edges, so a score computed with a different flank is not comparable to a published one.

**`--threads`** — The number of worker processes. Parallelism is at the chromosome level, so setting `--threads` above the number of chromosomes provides no benefit, and the largest chromosome sets the floor on run time. Each worker holds one chromosome's reference sequence, so memory grows with the worker count. The separate TSS pass batches overlapping windows (up to approximately 1 Mb of centre span) to reduce repeated BAM decompression. Each original TSS still contributes independently, including repeated annotations. The mate cache evicts unmatched records once their declared mate coordinate has passed; `diagnostics.pending_mates_peak` records the maximum per-chromosome cache size.


## Schema 2.0: definitions and missing values

Schema 2.0 separates unrestricted global editing counts from context statistics.
Global and per-fragment editing use all merged, quality-passing A/C/G/T calls at
reference C/G, including contig edges and sites adjacent to ambiguous reference
bases. `context_summary` counts only sites with a complete ACGT trinucleotide;
its totals need not equal `editing` totals. Older reports used context-restricted
global counts and should not be pooled with schema 2.0 without recomputation.

Means and medians of edit counts and fragment lengths use exact, untruncated
frequency tables. Display histograms can still have overflow bins. The
per-fragment rate median remains approximate (0.005-wide bins), with its method
recorded in JSON. Undefined rates and summaries are `null`, not zero; a zero TSS
background produces a null score and empty normalized CSV cells. JSON never
contains NaN or Infinity. `status` records unavailable calculations and TSS
annotation exclusions (unknown contigs and windows outside contig bounds).
Point TSS records with `start == end` are supported without shifting coordinates,
as are 1-bp and wider intervals. Negative starts, reversed intervals and malformed
coordinates fail before the expensive BAM scan; parsed sites are reused.
BAM/FASTA contig-length mismatches also fail early.

`file_reads` contains full-file indexed mapped/unmapped/total counts, including
unplaced unmapped reads. `reads` contains sampled record statistics; filtering
reason counts overlap and must not be added to derive a discarded total.
`provenance` records version, sample, input paths, timestamp and parameters.

## Deaminase and alignment diagnostics

- `editing.zero_edit_fraction` includes all analyzed fragments, including those
  without editable bases. The `zero_edit_by_opportunities.csv` stratifies by
  **0**, **1–10**, **11–25**, **26–50**, **51–100**, and **101+** opportunities.
  A zero-edit fragment is not by itself evidence of closed chromatin.
- `edit_directions.csv` is a sparse joint frequency table with columns
  `ct_edits,ga_edits,fragments`. Both directions are independent of read strand.
  JSON also reports direction totals and fragments containing both patterns.
- `diagnostics` reports non-edit mismatch frequency among merged, quality-filtered
  ACGT aligned bases, plus overlap disagreements and equal-quality conflicts
  relative to overlapping quality-passing positions. Equal-quality conflicts are
  removed from base counting.
- CIGAR insertion and soft-clip rates use passing-read query bases (M/I/S/=/X);
  deletion rate uses reference M/D/=/X bases. These are record-level measurements,
  include mate overlaps, and are not base-quality filtered.

`motif_enrichment.csv` records edited and opportunity counts/frequencies and their
ratio for each offset/base. Both sets use complete ACGT windows, the same filters,
and a common C-centred orientation. No pseudocount is added: missing background
or no edited events gives an empty ratio. This auxiliary CSV remains available
for downstream analysis, but its heatmap is not displayed in the HTML report.
The report retains the frequency/bits motif logo and trinucleotide context-rate
plot. These describe sequence preference; they do not perform footprint bias
correction. `context.csv` and `fragment_length.csv`
provide the numbers behind the corresponding plots. All diagnostic CSVs are
written with `--no_plot` as well.

## Compare multiple samples

```bash
uv run deamtools qc-summary --reports results/sample1.json results/sample2.json \
  --out_dir results/comparison --out_name cohort
```

Writes `cohort.csv` and a standalone `cohort.html` with sample identity, schema,
filters, sampling fraction, editing metrics, TSS, mismatch rates and context
preferences. Missing fields remain empty/n/a. Legacy reports are labelled;
no values are pooled and no automatic quality thresholds are applied.


## Interpretation rules

**Data checks PASS/FAIL refers only to data availability and arithmetic
consistency**, not to a biological assay-quality grade. The four named checks
require passing reads, editable bases, passing ≤ total reads, and
0 ≤ edited ≤ editable bases. An unavailable or inconsistent input fails these
checks; annotation exclusions, an unavailable requested TSS score, absent edits,
or observed variation across context rates generate advisory notes. The report
lists every check and advisory so the summary can be audited. Advisory notes
are displayed separately and never turn passing data checks into WARNING.
Green is reserved for validation, teal for descriptive observations, amber for
advisories, red for failed validation, and gray for unavailable values.

No universal edit-rate, duplicate-rate or TSS “ideal” threshold is configured.
Cards therefore use descriptive labels such as “Assay-dependent”, “Available”,
“Duplicate-flagged records” and “flank-normalized peak”. Zero duplicate flags
may reflect upstream deduplication and do not establish original library
complexity. TSS has a **1× flank-background reference**, not an invented acceptable/ideal gauge. This implementation's
unsmoothed 1-bp cut profile is not numerically identical to ENCODE's smoothed
coverage calculation.

The **Top-context / pooled edit rate** divides its edit fraction by the opportunity-weighted
pooled fraction of valid trinucleotide contexts. Both populations use the same
context constraints; the global editing rate is also displayed separately.
This is descriptive, not a statistical test, and no entropy or KL score is
introduced. Inspect opportunities for sparse contexts before interpretation.
The read-retention module also shows available filtering reason counts, including
QC-fail, low MAPQ, secondary and duplicate flags. Filter categories overlap and
do not sum to excluded reads; the dashboard does not invent a residual “other”
category or subtract overlapping flags sequentially.

## Portable directory output

```bash
uv run deamtools qc --bam sample.bam --fasta genome.fa \
  --out_dir results --out_name sample --report-dir results/sample_dashboard
```

Open `results/sample_dashboard/report.html`. Share the whole directory, including
`assets/`, CSVs and `metrics.json`. Large PNG and CSV payloads are external;
compact interactive SVG markup remains inline and an exportable copy is also
saved under `assets/`. Without `--report-dir`, PNGs and CSV downloads are embedded
in the single HTML. `--no_plot` retains all tables, status explanations and CSVs.

To rebuild presentation from saved outputs without rescanning BAM:

```python
import json
from deamtools.qc.report import write_report

with open("results/sample.json") as stream:
    metrics = json.load(stream)
write_report(metrics, "results", "sample",
             report_dir="results/sample_dashboard")
```

Only recorded fields and available companion CSVs are used. A legacy report
without a fragment-length CSV cannot reconstruct its original length plot;
summary values remain available. Rebuilding the report does not update the
underlying QC metrics or their schema.


### Shareable reports

```bash
uv run deamtools qc --bam sample.bam --fasta genome.fa \
  --out_dir results --out_name sample --redact-paths \
  --report-dir results/sample_shareable
```

The HTML source contains basenames for BAM, FASTA and TSS annotation, including
inside tooltips and collapsed sections. Directory components in the recorded
command and path-valued run parameters are removed as well; an executable such
as `/Users/name/.venv/bin/deamtools` becomes `deamtools`. Unparseable recorded
commands are omitted rather than copied verbatim. The portable `metrics.json`
is also redacted. The original `results/sample.json` and returned metrics retain
full provenance for local reproducibility. Sample names and basenames are retained;
this is path redaction, not complete anonymization.

When regenerating from saved outputs, pass `redact_paths=True` to `write_report`.
The input metrics dictionary is not modified. Default reports retain full paths.
