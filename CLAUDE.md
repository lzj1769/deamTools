# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DeamTools is a Python bioinformatics toolkit for deaminase-based chromatin accessibility analysis. It detects and quantifies cytosine-to-thymine (C→T) deamination events in aligned genomic reads, then converts them into BigWig coverage tracks for downstream epigenomics analysis.

**Status**: Early alpha (v0.1.0) — core `bam2bw` pipeline is partially implemented; several modules are stubs awaiting implementation.

## Setup

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.10. Key runtime dependencies: `numpy`, `pandas`, `pysam`, `pyBigWig`, `MOODS-python`. Optional: `pyjaspar` (for JASPAR motif fetching).

## Common Commands

```bash
# CLI entry point
deamtools --help
deamtools bam2bw --bam sample.bam --fasta hg38.fa --chrom_sizes hg38.sizes --output sample.bw

# Tests (no test files exist yet)
pytest
pytest --cov=deamtools

# Code quality
ruff check src/
black src/
mypy src/
```

## Architecture

Source lives under `src/deamtools/` (importable as `deamtools`). The CLI (`cli/main.py`) uses argparse and dispatches to subcommand implementations.

### Module Layout

| Module | Purpose |
|---|---|
| `cli/main.py` | Argparse dispatcher; currently wires `bam2bw` subcommand |
| `preprocessing/bam2bw.py` | BAM → BigWig conversion — main algorithm (incomplete, raises `NotImplementedError`) |
| `preprocessing/bam2fragment.py` | Fragment extraction stub |
| `preprocessing/fragment2bw.py` | Fragment → BigWig stub |
| `motif/_matching.py` | MOODS-based motif scanning; `prepare_scanner()` and JASPAR fetch are functional |
| `utils/chromosome.py` | Extract chrom sizes from BAM header or UCSC `.sizes` files |
| `utils/_logging.py` | Timestamped logging setup |
| `utils/version.py` | Dynamic version via `importlib.metadata` |

### Data Flow (intended)

```
BAM (coordinate-sorted) + FASTA reference
        ↓  bam2bw
  per-base C→T deamination counts
        ↓
  BigWig output (strand-aware coverage)
```

The `--regions` BED argument restricts processing to specific genomic intervals; `--extend_size` symmetrically extends each called editing site before writing signal.

## Known Issues / Development Notes

- `preprocessing/bam2bw.py`: references `pyBigWig` without importing it; variable `grs` is undefined on line 41.
- `motif/_matching.py`: contains a copy-pasted `bam2bw`-related block that does not belong there.
- No test files exist yet despite full pytest configuration.
