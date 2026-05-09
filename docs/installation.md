# Installation

## Requirements

| Requirement | Version |
|---|---|
| Python | ≥ 3.10 |
| samtools | any recent version |

samtools is required to sort and index BAM files and to index FASTA files before running DeamTools. Install it via your package manager or from [htslib.org](https://www.htslib.org/).

```bash
# macOS
brew install samtools

# Ubuntu / Debian
sudo apt-get install samtools

# conda
conda install -c bioconda samtools
```

## Install DeamTools

### From GitHub (recommended)

```bash
git clone https://github.com/lzj1769/deamTools.git
cd deamTools
pip install .
```

### Development install

If you plan to modify the source or run the test suite, install in editable mode with the `dev` extras:

```bash
pip install -e ".[dev]"
```

This additionally installs: `pytest`, `pytest-cov`, `ruff`, `black`, `mypy`.

### Documentation extras

To build the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs serve   # live-preview at http://127.0.0.1:8000
```

## Python dependencies

These are installed automatically by pip:

| Package | Purpose |
|---|---|
| `numpy` | Per-base count arrays and signal convolution |
| `pandas` | Tabular data utilities |
| `matplotlib` | Plotting utilities |
| `pysam` | BAM and FASTA file I/O |
| `pyBigWig` | BigWig file writing |
| `MOODS-python` | Motif scanning (future functionality) |

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

All 25 tests should pass. The test suite uses synthetic BAM and FASTA files created in a temporary directory — no external data files are required.
