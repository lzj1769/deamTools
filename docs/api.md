# API Reference

DeamTools can be used as a Python library as well as a command-line tool. All public functions are importable from their respective submodules.

## `deamtools.preprocessing.bam2bw`

### `run_bam2bw`

```python
from deamtools.preprocessing.bam2bw import run_bam2bw

run_bam2bw(
    bam_path: str,
    fasta_path: str,
    output_path: str,
    chrom_sizes_path: str | None = None,
    bed_path: str | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    extend_size: int = 0,
    threads: int = 1,
) -> None
```

Convert a BAM file to a BigWig track of per-base C→T deamination counts.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bam_path` | `str` | — | Path to a coordinate-sorted, indexed BAM file. |
| `fasta_path` | `str` | — | Path to the reference FASTA file (must have a `.fai` index). |
| `output_path` | `str` | — | Destination path for the output BigWig file. Parent directories are created automatically. |
| `chrom_sizes_path` | `str \| None` | `None` | Path to a tab-delimited chromosome sizes file. When `None`, sizes are read from the BAM header. |
| `bed_path` | `str \| None` | `None` | Path to a BED file of regions to restrict processing to. When `None`, the whole genome is processed. |
| `min_mapq` | `int` | `20` | Minimum read mapping quality. Reads below this value are skipped. |
| `min_baseq` | `int` | `20` | Minimum base quality at a position. Bases below this value are not counted. |
| `extend_size` | `int` | `0` | Symmetrically extend each deamination event by this many base pairs using box-kernel convolution. |
| `threads` | `int` | `1` | Number of threads for parallel chromosome processing. |

**Returns:** `None`. Writes a BigWig file to `output_path`.

**Example**

```python
from deamtools.preprocessing.bam2bw import run_bam2bw

run_bam2bw(
    bam_path="sample.sorted.bam",
    fasta_path="hg38.fa",
    output_path="results/sample.bw",
    min_mapq=30,
    min_baseq=30,
    threads=4,
)
```

---

### `_count_deamination_on_chrom`

```python
from deamtools.preprocessing.bam2bw import _count_deamination_on_chrom

chrom, counts = _count_deamination_on_chrom(
    bam_path: str,
    fasta_path: str,
    chrom: str,
    chrom_size: int,
    regions: list[tuple[int, int]] | None,
    min_mapq: int,
    min_baseq: int,
    extend_size: int,
) -> tuple[str, numpy.ndarray]
```

Low-level per-chromosome counting kernel. Opens its own BAM and FASTA handles (thread-safe) and returns a float32 array of deamination counts indexed by reference position.

Useful for custom workflows where you want per-chromosome count arrays without writing a BigWig.

**Parameters**

| Parameter | Type | Description |
|---|---|---|
| `bam_path` | `str` | Path to the indexed BAM file. |
| `fasta_path` | `str` | Path to the indexed FASTA file. |
| `chrom` | `str` | Chromosome name (must match the BAM and FASTA). |
| `chrom_size` | `int` | Length of the chromosome in base pairs. |
| `regions` | `list[tuple[int,int]] \| None` | List of `(start, end)` intervals to restrict counting to, or `None` for the whole chromosome. Intervals should be non-overlapping (use `_load_regions` to merge). |
| `min_mapq` | `int` | Minimum read mapping quality. |
| `min_baseq` | `int` | Minimum base quality at a position. |
| `extend_size` | `int` | Signal extension in base pairs (0 = no extension). |

**Returns:** `(chrom, counts)` where `counts` is a `numpy.ndarray` of shape `(chrom_size,)` and dtype `float32`.

**Example**

```python
import pysam
from deamtools.preprocessing.bam2bw import _count_deamination_on_chrom

with pysam.AlignmentFile("sample.bam", "rb") as bam:
    chrom_size = dict(zip(bam.references, bam.lengths))["chr1"]

chrom, counts = _count_deamination_on_chrom(
    bam_path="sample.sorted.bam",
    fasta_path="hg38.fa",
    chrom="chr1",
    chrom_size=chrom_size,
    regions=[(1_000_000, 2_000_000)],
    min_mapq=20,
    min_baseq=20,
    extend_size=0,
)

print(f"Total events on {chrom}: {int(counts.sum())}")
print(f"Events in first 1 Mb: {int(counts[:1_000_000].sum())}")
```

---

### `_load_regions`

```python
from deamtools.preprocessing.bam2bw import _load_regions

regions = _load_regions(bed_path: str) -> dict[str, list[tuple[int, int]]]
```

Parse a BED file and return a dictionary mapping chromosome names to sorted, non-overlapping intervals. Overlapping intervals on the same chromosome are merged.

**Parameters**

| Parameter | Type | Description |
|---|---|---|
| `bed_path` | `str` | Path to a BED file (tab-delimited; columns: chrom, start, end). Comment lines starting with `#` are skipped. |

**Returns:** `dict[str, list[tuple[int, int]]]` — keys are chromosome names; values are sorted lists of `(start, end)` tuples with no overlaps.

---

## `deamtools.utils`

### `get_chrom_sizes_from_bam`

```python
from deamtools.utils import get_chrom_sizes_from_bam

sizes = get_chrom_sizes_from_bam(bam: pysam.AlignmentFile) -> dict[str, int]
```

Extract chromosome names and lengths from an open BAM file's header.

**Parameters:** `bam` — an open `pysam.AlignmentFile` object.

**Returns:** `dict[str, int]` mapping chromosome names to their lengths.

---

### `get_chrom_sizes_from_file`

```python
from deamtools.utils import get_chrom_sizes_from_file

sizes = get_chrom_sizes_from_file(chrom_size_file: str) -> dict[str, int]
```

Parse a UCSC-style tab-delimited chromosome sizes file.

**Parameters:** `chrom_size_file` — path to a file where each line is `<chrom>\t<size>`.

**Returns:** `dict[str, int]` mapping chromosome names to their lengths.

---

### `get_version`

```python
from deamtools.utils import get_version

version = get_version() -> str
```

Return the installed package version string (e.g. `"0.1.0"`), read from the package metadata via `importlib.metadata`.
