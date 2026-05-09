# DeamTools

**DeamTools** is a Python command-line toolkit for deaminase-based chromatin accessibility analysis. It converts aligned sequencing reads into per-base BigWig coverage tracks of cytosine deamination events, enabling genome-wide quantification of chromatin accessibility at single-base resolution.

## What is deaminase-based chromatin accessibility?

Deaminase-based assays exploit the fact that deaminase enzymes preferentially act on single-stranded DNA, which is exposed in accessible (nucleosome-free) chromatin regions. The enzyme converts cytosine (C) to uracil (U), which is subsequently read as thymine (T) by the sequencer. After alignment, positions where the sequencing read carries T at a reference C site report an accessible region.

DeamTools automates the detection and quantification of these C→T editing events from BAM files, producing BigWig tracks ready for genome-browser visualisation, peak calling, or footprinting analysis.

## Key features

- **Single-base resolution** — counts deamination events at individual cytosine positions, not binned windows
- **Strand-aware** — detects C→T on forward-strand reads and G→A on reverse-strand reads (template-strand deamination)
- **Quality filtering** — configurable MAPQ and per-base quality thresholds; secondary, duplicate, QC-failed, supplementary reads are always excluded
- **Region-restricted mode** — process only BED-defined intervals instead of the whole genome, substantially reducing runtime for targeted analyses
- **Signal extension** — optionally extend each deamination site symmetrically to produce smoothed or windowed accessibility tracks
- **Multi-threaded** — chromosomes are processed in parallel for faster whole-genome runs

## Overview

```
Aligned reads (BAM)  +  Reference genome (FASTA)
            │
            ▼
   DeamTools bam2bw
   ┌─────────────────────────────────────┐
   │  For each read:                     │
   │    • skip low-quality / flagged     │
   │    • scan aligned base pairs        │
   │    • record C→T (fwd) / G→A (rev)  │
   │  Accumulate per-base counts         │
   │  Optionally extend signal           │
   │  Write sparse BigWig                │
   └─────────────────────────────────────┘
            │
            ▼
     Accessibility track (BigWig)
```

## Getting started

See the [Installation](installation.md) page to set up DeamTools, then the [bam2bw](usage/bam2bw.md) page for a full usage guide.
