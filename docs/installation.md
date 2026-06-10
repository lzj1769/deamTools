# Installation

## Requirements

| Requirement | Version | Needed for |
|---|---|---|
| Python | ≥ 3.10 | everything |
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

### From GitHub (recommended)

```bash
git clone https://github.com/lzj1769/deamTools.git
cd deamTools
pip install .
```

### Development install

To modify the source or run the test suite, install in editable mode with the `dev` extras:

```bash
pip install -e ".[dev]"
```

This additionally installs `pytest`, `pytest-cov`, `ruff`, `black`, and `mypy`.

### Documentation extras

This documentation is built with Sphinx (MyST Markdown + the Read the Docs theme). To build it locally:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html
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
| `MOODS-python` | Motif scanning |

## Verify installation

```bash
deamtools --version
# deamtools 0.1.0

deamtools --help
```

## Running tests

```bash
cd deamTools
pytest
```

The test suite uses synthetic BAM and FASTA fixtures created in a temporary directory, so no external data files — and no `bwa`/`samtools` — are required to run it.
