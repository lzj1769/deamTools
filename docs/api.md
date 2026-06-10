# API reference

Every subcommand is a thin wrapper around a public `run_*` function, so DeamTools can be used as a Python library as well as from the command line. Each function is importable from its submodule.

## `deamtools.align.index`

### `run_index`

```python
from deamtools.align.index import run_index

run_index(
    fasta_path: str,
    out_dir: str | None = None,
    out_name: str | None = None,
    force: bool = False,
) -> None
```

Build the deamination-aware reference index: the standard `<fasta>.fai` (always next to the FASTA) plus the doubly-converted reference and its BWA index at `<out_dir>/<out_name>.deamtools.c2t*`. `out_dir`/`out_name` default to the FASTA's own directory and file name. Requires `bwa` and `samtools`.

## `deamtools.align.align`

### `run_align`

```python
from deamtools.align.align import run_align

run_align(
    fasta_path: str,
    fastq1: str,
    out_dir: str,
    out_name: str,
    fastq2: str | None = None,
    threads: int = 1,
    read_group: str | None = None,
    index_path: str | None = None,
) -> None
```

Align deaminated reads (read 1 C→T, read 2 G→A) to a `deamtools index`-prepared reference and write a sorted, indexed `<out_dir>/<out_name>.bam`. Pass `index_path` when the index was built with a custom `--out_dir`/`--out_name`. Requires `bwa` and `samtools`.

## `deamtools.preprocessing.bam2bw`

### `run_bam2bw`

```python
from deamtools.preprocessing.bam2bw import run_bam2bw

run_bam2bw(
    bam_path: str,
    fasta_path: str,
    out_dir: str,
    out_name: str,
    chrom_sizes_path: str | None = None,
    bed_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    extend_size: int = 0,
    threads: int = 1,
    mode: str = "count",
    min_coverage: int = 10,
) -> None
```

Write a per-base BigWig of deamination signal to `<out_dir>/<out_name>.bw`. `mode="count"` writes raw edit counts (strand-agnostic `C→T` or `G→A`) and honours `extend_size`; `mode="ratio"` writes `edits / total_ACGT_coverage`, masking positions below `min_coverage` to 0.

## `deamtools.preprocessing.bam2fragment`

### `run_bam2fragment`

```python
from deamtools.preprocessing.bam2fragment import run_bam2fragment

run_bam2fragment(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
    threads: int = 1,
    barcode: bool = False,
    barcode_tag: str = "CB",
) -> None
```

Write a per-fragment editing-signal table to `output_path` (gzip-compressed when the path ends in `.gz`). With `barcode=True`, a barcode column is added (10x ordering). Editing is strand-aware (`C→T` forward, `G→A` reverse).

## `deamtools.qc`

### `run_qc`

```python
from deamtools.qc import run_qc

run_qc(
    bam_path: str,
    fasta_path: str,
    out_dir: str,
    out_name: str,
    tss_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    threads: int = 1,
    tss_flank: int = 2000,
    plot: bool = True,
) -> dict
```

Compute QC metrics and write `<out_dir>/<out_name>.json` plus a self-contained HTML report. Returns the metrics dictionary. Supplying `tss_path` adds TSS enrichment.

## `deamtools.stat.plot_motif`

### `run_plot_motif`

```python
from deamtools.stat.plot_motif import run_plot_motif

run_plot_motif(
    bigwig_path: str,
    fasta_path: str,
    output_path: str,
    bed_path: str | None = None,
    window_size: int = 10,
) -> pandas.DataFrame
```

Build a deaminase flanking-preference sequence logo from a count-mode BigWig, saving the plot to `output_path` and a sibling `.csv` of the bit-score matrix. Returns the bit-score matrix.

## `deamtools.utils`

```python
from deamtools.utils import (
    get_chrom_sizes_from_bam,   # (bam: pysam.AlignmentFile) -> dict[str, int]
    get_chrom_sizes_from_file,  # (path: str) -> dict[str, int]
    get_version,                # () -> str
)
from deamtools.utils.regions import _load_regions  # (bed_path: str) -> pandas.DataFrame
```

| Function | Returns |
|---|---|
| `get_chrom_sizes_from_bam(bam)` | Chromosome name → length, from an open BAM header. |
| `get_chrom_sizes_from_file(path)` | Chromosome name → length, from a UCSC `chrom.sizes` file. |
| `get_version()` | Installed package version (via `importlib.metadata`). |
| `_load_regions(bed_path)` | BED as a `chrom/start/end` DataFrame, with overlapping intervals merged. |
