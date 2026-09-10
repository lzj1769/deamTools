# Installation

## Requirements

| Requirement | Version | Needed for |
|---|---|---|
| Python | ≥ 3.12 | everything |
| samtools | any recent | indexing FASTA/BAM; used by `index` and `align` |
| bwa | any recent | `index` and `align` only |

`samtools` and `bwa` must be on your `PATH` for the `index` and `align` commands. The signal/QC commands (`bam2bw`, `bam2fragment`, `qc`) do not require them.

```bash
# macOS (Homebrew)
brew install samtools bwa

# Ubuntu / Debian
sudo apt-get install samtools bwa

# conda
conda install -c bioconda samtools bwa
```

## Install DeamTools

### From PyPI (recommended)

```bash
pip install deamtools
# or
uv pip install deamtools
```

To also install the optional `seq2edit` model, which pulls in PyTorch:

```bash
pip install "deamtools[seq2edit]"
```

### From GitHub

```bash
pip install git+https://github.com/lzj1769/deamTools.git
```

### Development install

The project is managed with [uv](https://docs.astral.sh/uv/), and `uv.lock` is
committed, so `uv sync` reproduces the development environment exactly:

```bash
git clone https://github.com/lzj1769/deamTools.git
cd deamTools
uv sync --extra dev
uv run pytest
```

This additionally installs `pytest`, `pytest-cov`, `ruff`, `black`, and `mypy`. Run
tools through `uv run` (for example `uv run ruff check src/ tests/`); there is no
need to activate the environment.

### Documentation extras

This documentation is built with Sphinx (MyST Markdown + the Read the Docs theme). To build it locally:

```bash
uv sync --extra docs
uv run sphinx-build -b html docs docs/_build/html
# open docs/_build/html/index.html
```

## Python dependencies

Installed automatically by pip:

| Package | Purpose |
|---|---|
| `numpy` | Per-base count arrays and signal convolution |
| `pandas` | Tabular / BED data utilities |
| `matplotlib` | Plotting (QC report, motif logo) |
| `logomaker` | Deaminase motif logo in the `qc` report |
| `pysam` | BAM and FASTA I/O |
| `pyBigWig` | BigWig reading and writing |
| `MOODS-python` | Motif scanning (`match`) |

The `match` command additionally needs the optional `pyjaspar` package to fetch motifs from JASPAR: `pip install pyjaspar`.

## Verify installation

```bash
deamtools --version
# deamtools 0.1.1

deamtools --help
```

## Running tests

```bash
cd deamTools
uv run pytest
```

The test suite uses synthetic BAM and FASTA fixtures created in a temporary directory, so no external data files — and no `bwa`/`samtools` — are required to run it.
