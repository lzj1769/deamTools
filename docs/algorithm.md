# Algorithm

This page describes how DeamTools detects and quantifies deamination events from aligned reads.

## Biological basis

Deaminase enzymes convert cytosine (C) to uracil (U) in single-stranded DNA. In deaminase-based chromatin accessibility assays, the enzyme is applied to chromatin in situ. Nucleosome-free (accessible) regions expose single-stranded DNA, making them preferred substrates for the deaminase. After library preparation and sequencing, deaminated cytosines appear as C→T substitutions in the reads relative to the reference genome.

The key signal is therefore a mismatch: the reference carries C, but the aligned read carries T at the same position. The density of such mismatches across the genome reports chromatin accessibility.

## Strand convention

Reads can originate from either strand of the double-stranded genome.

**Forward-strand reads** (`is_reverse = False`):

```
Reference (+ strand):  5'—...C...—3'
                               ↓ deamination
Read sequence:         5'—...T...—3'
```

Detection criterion: `ref_base == 'C'` and `read_base == 'T'`.

**Reverse-strand reads** (`is_reverse = True`):

The deaminase acts on the minus (template) strand, where a reference G on the plus strand corresponds to a C on the minus strand. After deamination (C→U/T on minus strand) and PCR amplification, the plus strand at that position reads A instead of G. The BAM stores the reverse-strand read as the reverse complement of the sequenced bases, so the mismatch appears as G→A.

```
Reference (+ strand):  5'—...G...—3'
                                   ↑ (template strand has C here)
Template (− strand):   3'—...C...—5'
                               ↓ deamination
Template after:        3'—...T...—5'
After PCR / rev-comp:  5'—...A...—3'  ← stored in BAM
```

Detection criterion: `ref_base == 'G'` and `read_base == 'A'`.

Both strand types are detected and their events are accumulated into the same per-base count array indexed by reference position.

## Processing pipeline

### 1. Determine chromosomes to process

Chromosome names and sizes are taken from either a user-supplied sizes file or the BAM header. If a BED file of regions is provided, only chromosomes that appear in the BED file are processed.

### 2. Load and merge BED regions (if provided)

BED intervals are loaded per chromosome and sorted by start position. Overlapping intervals are merged into non-overlapping ones to prevent a read spanning an overlap from being counted twice:

```
Input:   [0, 500)  [300, 800)  [1000, 1500)
Merged:  [0, 800)              [1000, 1500)
```

### 3. Count deamination events per chromosome

For each chromosome, a float32 array of length equal to the chromosome size is initialised to zero. For each genomic region (the full chromosome, or each merged BED interval), the reference sequence is fetched once with `pysam.FastaFile`. The BAM is then iterated with `pysam.AlignmentFile.fetch()`, which uses the BAM index to skip regions with no reads.

For each read that passes all filters:

```
for (query_pos, ref_pos) in read.get_aligned_pairs(matches_only=True):
    if ref_pos outside region: skip
    if base quality < min_baseq: skip

    ref_base = reference[ref_pos]
    read_base = read_sequence[query_pos]

    if forward read and ref_base == C and read_base == T:
        counts[ref_pos] += 1
    if reverse read and ref_base == G and read_base == A:
        counts[ref_pos] += 1
```

`get_aligned_pairs(matches_only=True)` returns only positions that are aligned matches (no insertions, deletions, or soft clips), so the query and reference bases are always directly comparable.

### 4. Signal extension (optional)

When `--extend_size > 0`, the raw count array is convolved with a box kernel of width `2 × extend_size + 1`:

```
kernel = [1, 1, 1, ..., 1]   # 2 * extend_size + 1 ones
counts  = convolve(counts, kernel, mode='same')
```

This means the value at position *p* after convolution equals the number of raw deamination events in the window [*p* − extend_size, *p* + extend_size]. The convolution is computed with NumPy using zero-padding at chromosome boundaries.

### 5. Parallel execution

Each chromosome is processed in a separate thread via `concurrent.futures.ThreadPoolExecutor`. The BAM and FASTA files are opened independently within each thread, so there is no shared mutable state. Results (one array per chromosome) are collected after all threads finish, then written to the BigWig sequentially.

### 6. Write BigWig

Only positions with non-zero signal are written (sparse format), using `pyBigWig.addEntries` with `span=1`. Chromosomes are written in the same order they appear in the chromosome sizes dictionary (preserving BAM header order). The BigWig header lists all chromosomes, even those with no events.

## Read filtering

The following filters are applied before any base-level inspection:

| Filter | Condition |
|---|---|
| Unmapped | `read.is_unmapped` |
| PCR/optical duplicate | `read.is_duplicate` |
| QC fail | `read.is_qcfail` |
| Secondary alignment | `read.is_secondary` |
| Supplementary alignment | `read.is_supplementary` |
| Low mapping quality | `read.mapping_quality < min_mapq` |
| No sequence | `read.query_sequence is None` |

At the per-base level, positions where the base quality is below `min_baseq` are skipped individually; the rest of the read is still processed.

## Complexity

Let *N* be the number of aligned reads and *L* the average read length. The algorithm runs in O(*N* × *L*) time for a single thread. Memory usage per chromosome is O(*C*) where *C* is the chromosome size (one float32 per base, ~1 GB for the largest human chromosomes). With *T* threads, up to *T* chromosomes are held in memory simultaneously.
