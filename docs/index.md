# Introduction

**DeamTools** is an open-source Python 3.10+ command-line toolkit for the analysis of
**deaminase-based chromatin accessibility** data. Double-stranded DNA cytosine
deaminases (e.g. DddA, DddSs/SsdAtox) preferentially edit cytosines in accessible,
protein-free chromatin. DeamTools turns the resulting single-base editing events into
genome-wide accessibility and transcription-factor footprinting signal.

It is developed for assays such as ACCESS-ATAC and related deaminase footprinting
methods, and provides a small set of composable commands that cover the full workflow
from raw reads to quantitative tracks and quality reports.

## How it works

A deaminase converts cytosine (C) to uracil (U) on exposed single-stranded DNA;
after PCR and sequencing this is read as a C→T substitution (or G→A on the opposite
strand). The density of these edits reports chromatin accessibility, while proteins
bound to DNA shield it from editing, leaving footprints. DeamTools detects these
edits from aligned reads and summarises them at single-base resolution.

## Components

DeamTools is organised as a set of subcommands, run as `deamtools <command>`:

| Command | Purpose |
|---|---|
| [`index`](usage/index.md) | Build the deamination-aware reference index (FASTA `.fai` + converted BWA index). |
| [`align`](usage/align.md) | Align deaminated reads (bwa-meth-style) to the indexed reference, producing a sorted BAM. |
| [`bam2bw`](usage/bam2bw.md) | Convert a BAM to a per-base BigWig of editing counts or conversion ratios. |
| [`bam2fragment`](usage/bam2fragment.md) | Convert a BAM to a per-fragment editing-signal table (bulk or single-cell). |
| [`qc`](usage/qc.md) | Quality-control metrics with a self-contained HTML report. |
| [`plot_motif`](usage/plot_motif.md) | Build a deaminase sequence-preference logo from an editing-count BigWig. |

## Workflow

```
FASTQ ─▶ index ─▶ align ─▶ BAM ─┬─▶ bam2bw ──────▶ BigWig ─▶ plot_motif
                                ├─▶ bam2fragment ─▶ fragment table
                                └─▶ qc ──────────▶ JSON + HTML report
```

## Supported file formats

| Format | Used by | Role |
|---|---|---|
| FASTA (`.fa` + `.fai`) | all | Reference genome. |
| FASTQ (`.fq[.gz]`) | `align` | Raw sequencing reads. |
| BAM (`.bam` + `.bai`) | most | Coordinate-sorted aligned reads. |
| BigWig (`.bw`) | `bam2bw`, `plot_motif` | Per-base editing signal. |
| BED | `bam2bw`, `qc`, `plot_motif` | Regions / TSS to restrict analysis to. |

## Getting started

1. Install DeamTools and `samtools`/`bwa` — see [Installation](installation.md).
2. Index your reference, then align reads — see [index](usage/index.md) and [align](usage/align.md).
3. Generate signal and QC — see [bam2bw](usage/bam2bw.md) and [qc](usage/qc.md).

To understand how editing events are detected and counted, see the
[Algorithm](algorithm.md) page; for the Python API, see the
[API reference](api.md).

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
usage/plot_motif
```

```{toctree}
:hidden:
:caption: Reference

algorithm
api
```
