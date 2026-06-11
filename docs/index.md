# Introduction

**DeamTools** is an open-source Python 3.10+ command-line toolkit for the analysis of
**deaminase-based chromatin accessibility** data. Double-stranded DNA cytosine
deaminases (e.g. DddA, DddSs/SsdAtox) preferentially edit cytosines in accessible,
protein-free chromatin. DeamTools turns the resulting single-base editing events into
genome-wide accessibility and transcription-factor footprinting signal.

It is developed for assays such as ACCESS-ATAC and related deaminase footprinting
methods, and provides a small set of composable commands that cover the full workflow
from raw reads to quantitative tracks and quality reports.

## Background

Chromatin accessibility marks the regulatory regions of the genome — promoters,
enhancers, and the binding sites of transcription factors (TFs). The established
assays, **DNase-seq** and **ATAC-seq**, read out accessibility from the *ends* of
enzymatically cut or transposed fragments, which yields at most two informative
positions per fragment and limits resolution to roughly the fragment length.

Deaminase-based assays instead use a double-stranded DNA cytosine deaminase to
**write many marks along each accessible molecule**. Every editable cytosine in
exposed DNA can be converted, so a single read carries dozens of accessibility
measurements rather than two. This gives several advantages:

- **Single-base resolution** — editing is recorded per cytosine, not per fragment end.
- **Footprints** — DNA bound by a TF or nucleosome is shielded from editing, leaving a
  local depletion that pinpoints the bound element.
- **Single-molecule / single-allele readout** — because the marks are encoded as
  sequence changes, each read reports the accessibility state of one DNA molecule, and
  is compatible with standard short-read sequencing and PCR.

## How it works

A deaminase converts cytosine (C) to uracil (U) on exposed single-stranded DNA;
after PCR and sequencing this is read as a **C→T** substitution on the top strand, or
a **G→A** substitution on reads from the bottom strand. The density of these edits
reports chromatin accessibility, while proteins bound to DNA shield it from editing,
leaving footprints.

DeamTools handles the two computational challenges this creates. First, heavily edited
reads align poorly to an unmodified genome, so `index`/`align` use a bwa-meth-style
three-letter alignment against a converted reference. Second, the edits must be
distinguished from the reference and tallied; `bam2bw`, `bam2fragment`, and `qc`
detect C→T / G→A mismatches at aligned positions and summarise them — as per-base
tracks, per-fragment tables, or QC metrics — at single-base resolution.

## Key features

- **Deamination-aware alignment** of heavily converted reads, with original sequences
  and chromosome names restored in a standard sorted BAM.
- **Per-base editing tracks** (BigWig) as raw counts or conversion ratios, with optional
  signal extension and region restriction.
- **Per-fragment tables** capturing single-molecule edit patterns, with cell-barcode
  support for single-cell data.
- **Quality control** in a self-contained HTML report: editing rate, per-read edit-rate
  distribution, enzyme context bias, fragment sizes, and TSS enrichment.
- **Enzyme-bias diagnostics**, including a deaminase sequence-preference logo built
  directly from the BAM and embedded in the QC report.
- **Multi-threaded** and **region-restricted** processing for whole-genome or targeted runs.

## Components

DeamTools is organised as a set of subcommands, run as `deamtools <command>`:

| Command | Purpose |
|---|---|
| [`index`](usage/index.md) | Build the deamination-aware reference index (FASTA `.fai` + converted BWA index). |
| [`align`](usage/align.md) | Align deaminated reads (bwa-meth-style) to the indexed reference, producing a sorted BAM. |
| [`bam2bw`](usage/bam2bw.md) | Convert a BAM to a per-base BigWig of editing counts or conversion ratios. |
| [`bam2fragment`](usage/bam2fragment.md) | Convert a BAM to a per-fragment editing-signal table (bulk or single-cell). |
| [`qc`](usage/qc.md) | Quality-control metrics — including the deaminase sequence-motif logo — in a self-contained HTML report. |
| [`matching`](usage/matching.md) | Scan regions for TF motif matches (MOODS) and write a BED of binding sites. |

## Workflow

```
FASTQ ─▶ index ─▶ align ─▶ BAM ─┬─▶ bam2bw ──────▶ BigWig track
                                ├─▶ bam2fragment ─▶ fragment table
                                └─▶ qc ──────────▶ JSON + HTML report (+ motif logo)
```

## Quick start

A typical end-to-end run, from a reference and FASTQs to tracks and a QC report:

```bash
# 1. One-time: index the reference (writes <fasta>.fai + the converted BWA index)
deamtools index --fasta genome.fa

# 2. Align deaminated reads -> results/sample.bam (+ .bai)
deamtools align --fasta genome.fa \
    --read1 sample_R1.fq.gz --read2 sample_R2.fq.gz \
    --out_dir results --out_name sample

# 3. Per-base editing track -> results/sample.bw
deamtools bam2bw --bam results/sample.bam --fasta genome.fa \
    --out_dir results --out_name sample

# 4. Quality-control report (with the deaminase motif logo)
#    -> results/sample.json + results/sample.html
deamtools qc --bam results/sample.bam --fasta genome.fa \
    --out_dir results --out_name sample
```

See each [command page](usage/index.md) for the full option list.

## Supported file formats

| Format | Used by | Role |
|---|---|---|
| FASTA (`.fa` + `.fai`) | all | Reference genome. |
| FASTQ (`.fq[.gz]`) | `align` | Raw sequencing reads. |
| BAM (`.bam` + `.bai`) | most | Coordinate-sorted aligned reads. |
| BigWig (`.bw`) | `bam2bw` | Per-base editing signal. |
| BED | `bam2bw`, `qc` | Regions / TSS to restrict analysis to. |

## Getting started

1. Install DeamTools and `samtools`/`bwa` — see [Installation](installation.md).
2. Index your reference, then align reads — see [index](usage/index.md) and [align](usage/align.md).
3. Generate signal and QC — see [bam2bw](usage/bam2bw.md) and [qc](usage/qc.md).

To understand how editing events are detected and counted, see the
[Algorithm](algorithm.md) page; for the Python API, see the
[API reference](api.md).

## Related methods and further reading

Deaminase-based accessibility and footprinting is an active area; DeamTools targets
ACCESS-ATAC but the file formats and metrics apply broadly. Key methods:

- **ACCESS-ATAC** — Yu, Li, *et al.* Deaminase-mediated chromatin accessibility profiling
  with single-allele resolution. *bioRxiv* (2024).
  [doi:10.1101/2024.12.17.628768](https://doi.org/10.1101/2024.12.17.628768)
- **DAF-seq** — Swanson *et al.* Mapping single-cell diploid chromatin fiber architectures.
  *Nature Biotechnology* (2025).
  [doi:10.1038/s41587-025-02914-3](https://doi.org/10.1038/s41587-025-02914-3)
- **TDAC-seq** — Roh *et al.* Coupling CRISPR scanning with targeted chromatin accessibility
  profiling using a double-stranded DNA deaminase. *Nature Methods* (2025).
  [doi:10.1038/s41592-025-02811-2](https://doi.org/10.1038/s41592-025-02811-2)
- **cFOOT-seq** — Wang *et al.* Genome-wide investigation of transcription factor footprints
  using cFOOT-seq. *Protein & Cell* (2025).
  [doi:10.1093/procel/pwaf071](https://doi.org/10.1093/procel/pwaf071)
- **FOODIE** — He *et al.* Genome-wide single-cell and single-molecule footprinting of
  transcription factors with deaminase. *PNAS* (2024).
  [doi:10.1073/pnas.2423270121](https://doi.org/10.1073/pnas.2423270121)

## Source and issues

DeamTools is developed on GitHub at
[github.com/lzj1769/deamTools](https://github.com/lzj1769/deamTools). Please report
bugs and feature requests via the issue tracker. Released under the MIT License.

```{toctree}
:hidden:
:caption: Getting started

installation
```

```{toctree}
:hidden:
:caption: Commands

usage/index
usage/align
usage/bam2bw
usage/bam2fragment
usage/qc
usage/matching
```

```{toctree}
:hidden:
:caption: Reference

algorithm
api
```
