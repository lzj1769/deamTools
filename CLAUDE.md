# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DeamTools is a Python command-line toolkit for deaminase-based chromatin accessibility analysis. In these assays a deaminase marks accessible chromatin by converting cytosine to uracil/thymine; the resulting single-base C→T (forward strand) / G→A (reverse strand) editing events are the accessibility signal. The toolkit covers the full path from raw reads to QC: deamination-aware alignment → per-base BigWig tracks / per-fragment tables → QC report (with a deaminase motif logo) → TF-motif matching.

The package root is `deamTools/` (note the capital T); the importable package is `deamtools` (lowercase) under `src/`. The outer directory also holds `data/` (sample BAMs/peaks/FASTAs, gitignored) and `test/` (ad-hoc shell scratch scripts, **not** the test suite).

## Commands

```bash
# Install (editable, with dev tooling) — run from deamTools/
pip install -e ".[dev]"

# Tests (pure-Python; addopts = -q, pythonpath = src)
pytest
pytest tests/test_qc.py                 # single file
pytest tests/test_qc.py::TestQC::test_edit_rate_and_opportunities  # single test

# Lint / format / types (config in pyproject.toml; line-length 88)
ruff check src/ tests/
black src/ tests/
mypy src/

# Docs (Sphinx + MyST + sphinx_rtd_theme)
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html      # CI uses -W (warnings = errors)
```

Tests synthesize BAM/FASTA/BigWig fixtures with `pysam`/`pyBigWig` in tmp dirs (see helpers like `_make_read`/`_write_bam` in `tests/test_bam2bw.py`). They do **not** require `bwa`/`samtools` on PATH; `MOODS` is a hard dependency and is exercised directly (`tests/test_matching.py` builds motifs in memory). The `index`/`align` commands need `bwa`+`samtools` at runtime; `match` needs the optional `pyjaspar` for JASPAR fetch.

## Architecture

CLI entry point is `deamtools.cli.main:main` (also `python -m deamtools`). `cli/main.py` is a single argparse dispatcher: each subcommand has an `_add_*_parser` builder and an `_run_*` thunk that unpacks the namespace and calls the module's public `run_*` function. `_log_invocation` echoes resolved arguments MACS2-style at INFO. **All CLI surface lives in `cli/main.py`**; implementation modules expose only the `run_*` entry point plus private helpers.

| Command | Module | What it does |
|---|---|---|
| `index` | `align/index.py` | `run_index` — builds the doubly-converted BWA index (+ FASTA `.fai`) |
| `align` | `align/align.py` | `run_align` — converts reads, runs `bwa mem`, restores, sorts to BAM |
| `bam2bw` | `preprocessing/bam2bw.py` | `run_bam2bw` — per-base editing **count** or **ratio** BigWig |
| `bam2fragment` | `preprocessing/bam2fragment.py` | `run_bam2fragment` — per-fragment editing-signal table (bulk or barcoded) |
| `qc` | `qc/qc.py` | `run_qc` — QC metrics → JSON + self-contained HTML report (incl. motif logo) |
| `match` | `motif/match.py` | `run_motif_matching` — MOODS motif scan → BED of binding sites |
| `footprint` | `footprint/footprint.py` | `run_footprint` — TF footprint score + permutation p-value at motif sites, from a BigWig |
| `seq2edit train` | `seq2edit/` | `run_train` — train a CNN (`EditNet`) mapping one-hot DNA → per-base editing rate (deaminase sequence bias), Poisson loss |
| — | `utils/` | `_load_regions` (BED parse+merge), `get_chrom_sizes_from_{bam,file}`, `get_version`, `logger` |

`preprocessing/fragment2bw.py` is an empty placeholder (not wired up).

### Output-path convention

All commands that write output (`index`, `align`, `bam2bw`, `bam2fragment`, `qc`, `match`) take `--out_dir` + `--out_name` and write `<out_dir>/<out_name>.<ext>` (`.bam`/`.bw`/`.tsv`/`.json`+`.html`/`.bed`; `bam2fragment --gzip` adds `.gz`). For `index`/`align` these are optional/derived where noted (`index` defaults the converted-index location next to the FASTA; `align --index` points at a custom one). The FASTA `.fai` is always written next to the FASTA (the pysam-based commands require it there).

### The deamination-aware alignment scheme (index + align)

The conceptual core, spanning two files (bwa-meth / three-letter strategy):

- **`index`** writes a *doubly-converted* reference (`<...>.deamtools.c2t`): every chromosome appears twice — a C→T copy prefixed `f` and a G→A copy prefixed `r` — then `bwa index`es it.
- **`align`** uses a **dual-conversion, take-best** strategy (the ACCESS-ATAC deaminase edits both strands, so one read can carry both `C→T` and `G→A`). `_feed_converted` emits each read twice — single-end as `ct`/`ga` candidates, paired-end as two fragment orientations `f` = (R1 C→T, R2 G→A) and `r` = (R1 G→A, R2 C→T) — sharing the read name, with the original SEQ in `YS:Z:` and the candidate in `YC:Z:` (carried through `bwa mem -C`). A feeder thread pumps the candidates into `bwa`'s stdin; the main thread groups bwa's output by read name (`_process_sam` relies on bwa preserving input order so a read's candidates are consecutive), keeps the higher-scoring candidate (`_primary_score` = sum of mates' `AS`), then `_restore_alignment` restores SEQ from `YS`, strips the `f`/`r` RNAME/RNEXT prefix, and drops `YS`/`YC`; `@HD`/`@SQ` are rebuilt from the FASTA `.fai`.

**Flow:** the restored SAM is written to `<out_dir>/<out_name>.sam`, then `samtools sort` converts it to a coordinate-sorted `<out_name>.bam` and it is indexed. (This is now a write-then-sort flow, not a single concurrent pipe; the `.sam` is left in place.)

### Signal computation conventions (shared across modules)

- **Edit definition.** `bam2bw` and `qc` count edits *strand-agnostically*: any reference `C→T` **or** `G→A` mismatch regardless of read orientation (see `_get_edit_count`; `qc._process_chrom`). `bam2fragment` is *strand-aware*: forward reads record `C→T`, reverse reads record `G→A`. Keep this distinction when changing either.
- **Read filtering.** Unmapped / duplicate / QC-fail / secondary / supplementary reads are always excluded, then `--min_mapq`; `--min_baseq` gates individual bases. Aligned positions come from `get_aligned_pairs(matches_only=True)`.
- **Parallelism.** `bam2bw`, `bam2fragment`, and `qc` fan out over regions/chromosomes with `ThreadPoolExecutor`; worker functions **open their own** `pysam` handles (pysam handles are not thread-safe).
- **BigWig ordering.** `pyBigWig` requires entries in `(chrom, start, end)` order; `run_bam2bw` sorts regions by BAM-header chromosome order and emits only non-zero bases (`span=1`).
- **`bam2bw` modes.** `count` writes raw edit counts and honors `--extend_size` (broadcasts each event into a `2*extend_size+1` window); `ratio` writes `edits / total_ACGT_coverage`, masks positions below `--min_coverage` to 0, and ignores `--extend_size`. `--normalize` (count mode only; warns and is ignored in `ratio`) scales every value by `--scale_factor / total` so the track sums to `--scale_factor` (default `1e6`); `total` is summed over the **processed regions after `--extend_size` broadcasting**, so with `--regions` it is not a genome-wide count. Because that total spans all regions, `run_bam2bw` holds every region's signal array in memory (the `results` dict) before writing — the one place memory scales with the region set rather than with a single region.
- **`qc`.** Single pass over the BAM accumulates: read stats, editing rate, per-read edit-rate distribution, trinucleotide context bias, fragment-length distribution, optional TSS enrichment, and a **deaminase motif PWM** (window around each edit, centre excluded, G-edits reverse-complemented). Outputs JSON + a self-contained HTML report that embeds a matplotlib summary figure and a `logomaker` motif logo (both base64, `Agg` backend).
- **`match`.** Builds a `MOODS.scan.Scanner` (`prepare_scanner`: log-odds matrices, p-value thresholds, reverse complements laid out as `[fwd_0..fwd_{n-1}, rc_0..rc_{n-1}]`), scans each region's sequence, and writes 6-column BED (`chrom start end motif score strand`). Motifs come from JASPAR (`pyjaspar`) or a caller-supplied list.
- **`footprint`.** For each motif site (width `L`) reads a per-base BigWig over `[start-L, end+L)` and scores `mean(left flank) + mean(right flank) - mean(centre)`; a positive score gets a permutation p-value from `n_shuffles` within-window shuffles (vectorised with a NumPy `Generator`). Parses the regions BED directly (not `_load_regions`, which merges and drops names) and parallelises over chromosomes. Writes `chrom start end name fp_score p_value`.
- **`seq2edit`.** The one **nested** subcommand group (`deamtools seq2edit <train|…>`, `required=True` subparser) and the one feature needing the **optional** `torch` extra (`pip install 'deamtools[seq2edit]'`). `model.py` holds both the torch-free encoding/data-prep helpers (`one_hot_encode`, window tiling, `build_dataset`, conv-size math `flatten_dim`) and `EditNet`; importing `model.py` pulls in torch (for `EditNet`), so `train.py` keeps its top-level imports torch-free and defers `EditNet`/`build_dataset`/`dataset` imports inside `run_train` (via `_require_torch`, which raises an actionable hint) — importing the CLI therefore never pulls in torch. `train` tiles `--train_regions` into `--seq_len` windows (input = one-hot ref sequence, target = per-base `--bigwig` signal, treated as Poisson counts), trains `EditNet` (2 conv blocks → FC head → per-base **Softplus** rate λ) with a **Poisson NLL** loss (`PoissonNLLLoss(log_input=False)`) + Adam + `ReduceLROnPlateau` (ACCESS-ATAC defaults), and writes the best-validation checkpoint `<out_name>.pth` (state dict + rebuild `config` + loss history). `EditNet` mirrors ACCESS-ATAC's `cnn_bias_model` but computes the FC input size from the architecture (`flatten_dim`) instead of hard-coding `928`, and uses a Softplus rate head + Poisson loss instead of the reference's ReLU + MSE. `model.py` imports torch, so `tests/test_seq2edit.py` gates the whole file with `pytest.importorskip("torch")`.

The counting/alignment logic intentionally mirrors the upstream ACCESS-ATAC reference (pinellolab/ACCESS-ATAC); preserve that correspondence.

## Conventions

- Python ≥ 3.10; every module uses `from __future__ import annotations`. Type hints expected but `disallow_untyped_defs = false`.
- A leading underscore marks a package-private symbol (e.g. `_load_regions`, `_get_edit_count`); some are imported across modules within the package, so treat them as package-internal rather than file-private.
- New public functions follow the existing NumPy-style docstring convention (see `bam2bw.py` / `qc.py`).
- `version.py` resolves the version via `importlib.metadata`; the canonical version lives in `pyproject.toml`.
- Docs are **Sphinx** (`docs/conf.py`, MyST Markdown pages, `sphinx_rtd_theme`); the GitHub Actions `docs` workflow builds them and publishes `docs/_build/html` to the `gh-pages` branch.

## Git / workflow

This repo (the `deamTools/` subdirectory) is the git root; remote `origin` is `github.com:lzj1769/deamTools`. Per the maintainer's preference, commit and push **directly to `main`** here (no feature branch / PR). End commit messages with a `Co-Authored-By: Claude <model>` trailer naming the model that made the change (currently Opus 5).
