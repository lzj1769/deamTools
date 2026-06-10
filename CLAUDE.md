# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DeamTools is a Python command-line toolkit for deaminase-based chromatin accessibility analysis. In these assays a deaminase marks accessible chromatin by converting cytosine to uracil/thymine; the resulting single-base C→T (forward strand) / G→A (reverse strand) editing events are the accessibility signal. The toolkit covers the full path from raw reads to visualization: deamination-aware alignment → per-base BigWig tracks or per-fragment tables → deaminase motif logos.

The package root is `deamTools/` (note the capital T); the importable package is `deamtools` (lowercase) under `src/`. The outer directory also holds `data/` (sample BAMs/peaks, gitignored, used for manual runs) and `test/` (ad-hoc shell scratch scripts, not the test suite).

## Commands

```bash
# Install (editable, with dev tooling) — run from deamTools/
pip install -e ".[dev]"

# Tests
pytest                              # all tests (addopts = -q)
pytest -v
pytest tests/test_bam2bw.py        # single file
pytest tests/test_bam2bw.py::test_name   # single test

# Lint / format / types (config in pyproject.toml; line-length 88)
ruff check src/ tests/
black src/ tests/
mypy src/
```

Tests are pure-Python: they synthesize BAM/FASTA/BigWig fixtures with `pysam`/`pyBigWig` in tmp dirs (see `tests/test_bam2bw.py` helpers `_make_read`/`_write_bam`). They do **not** require `bwa`/`samtools` on PATH — but the `index` and `align` subcommands do at runtime.

## Architecture

CLI entry point is `deamtools.cli.main:main` (also `python -m deamtools`). `cli/main.py` is a single argparse dispatcher: each subcommand has an `_add_*_parser` builder and an `_run_*` thunk that unpacks the namespace and calls the module's public `run_*` function. `_log_invocation` echoes the resolved arguments MACS2-style at INFO. **All CLI surface lives in `cli/main.py`**; the implementation modules expose only the `run_*` entry point plus private helpers.

Subcommands and their implementations:

| Command | Module | What it does |
|---|---|---|
| `index` | `align/index.py` | `run_index` — builds a doubly-converted BWA index next to the FASTA |
| `align` | `align/align.py` | `run_align` — converts reads on the fly, runs `bwa mem`, restores original sequences |
| `bam2bw` | `preprocessing/bam2bw.py` | `run_bam2bw` — per-base editing **count** or **ratio** BigWig |
| `bam2fragment` | `preprocessing/bam2fragment.py` | `run_bam2fragment` — per-fragment editing-signal table (bulk or barcoded) |
| `plot_motif` | `stat/plot_motif.py` | `run_plot_motif` — deaminase flanking-sequence logo from a count BigWig |
| — | `utils/` | `_load_regions` (BED parse+merge), `get_chrom_sizes_from_{bam,file}`, `get_version`, logging |

### The deamination-aware alignment scheme (index + align)

This is the conceptual core and spans two files. Following the bwa-meth strategy:

- **`index`** writes a *doubly-converted* reference (`<fasta>.deamtools.c2t`): every chromosome appears twice — a C→T copy prefixed `f` and a G→A copy prefixed `r` — then `bwa index`es it. The original `.fai` is also produced.
- **`align`** streams FASTQ through an on-the-fly converter (read 1 → C→T, read 2 → G→A), stashing the original SEQ in a `YS:Z:` tag carried through `bwa mem -C`. A post-processing pass over the SAM stream (`_restore_alignment`) restores the original SEQ from `YS`, strips the `f`/`r` RNAME/RNEXT prefix, and the `@HD`/`@SQ` header is regenerated from the original FASTA's `.fai`. Output is piped to `samtools sort` and indexed.

The pipeline runs as concurrent subprocesses (`bwa mem | post-process | samtools sort`) with a feeder thread pumping converted reads into `bwa`'s stdin; threads are split between bwa and sort.

### Signal computation conventions (shared across modules)

- **Edit definition.** `bam2bw` counts edits *strand-agnostically*: any reference `C→T` **or** `G→A` mismatch (see `_get_edit_count`). `bam2fragment` is *strand-aware*: forward reads record `C→T`, reverse reads record `G→A` (see `_editing_positions`). Keep this distinction in mind when changing either.
- **Read filtering.** Unmapped / duplicate / QC-fail / secondary / supplementary reads are always excluded, then `--min_mapq` is applied; `--min_baseq` gates individual bases. Aligned positions come from `get_aligned_pairs(matches_only=True)` so indels are handled naturally.
- **Parallelism.** Both preprocessing commands fan out over regions/chromosomes with `ThreadPoolExecutor`. Worker functions **open their own** `pysam` handles because pysam handles are not thread-safe.
- **BigWig ordering.** `pyBigWig` requires entries added in `(chrom, start, end)` order; `run_bam2bw` sorts regions by the BAM-header chromosome order before writing, and only non-zero bases are emitted (`span=1`).
- **`bam2bw` modes.** `count` writes raw edit counts and honors `--extend_size` (broadcasts each event into a `2*extend_size+1` window). `ratio` writes `edits / total_ACGT_coverage`, masking positions below `--min_coverage` to 0, and ignores `--extend_size`.
- **`plot_motif`.** Consumes a **count-mode, `--extend_size 0`** BigWig (a non-zero extension would skew the logo). It accumulates a PWM from windows around non-zero bases, excludes the center (editing) base so the logo reflects flanking preference, reverse-complements windows whose center is `G` to unify C→T/G→A, converts to bits, renders with `logomaker`, and writes a sibling `.csv`. Uses the `Agg` matplotlib backend (headless-safe).

The algorithm intentionally mirrors the upstream ACCESS-ATAC reference (pinellolab/ACCESS-ATAC); preserve that correspondence when modifying counting logic.

## Conventions

- Python ≥ 3.10; every module uses `from __future__ import annotations`. Type hints are expected but `disallow_untyped_defs = false`.
- A leading underscore marks a module-private function/symbol (e.g. `_load_regions`, `_get_edit_count`) — but note these are sometimes imported across modules within the package, so treat them as package-internal rather than file-private.
- New public functions follow the existing NumPy-style docstring convention (see `bam2bw.py` / `plot_motif.py`).
- `version.py` resolves the version dynamically via `importlib.metadata`; the canonical version lives in `pyproject.toml`.
- Docs are MkDocs (`docs/`, `mkdocs.yml`); install with `pip install -e ".[docs]"`.
