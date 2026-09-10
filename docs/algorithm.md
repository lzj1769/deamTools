# Algorithm

This page describes how DeamTools aligns deaminated reads and detects, quantifies, and summarises deamination events.

## Biological basis

Double-stranded DNA cytosine deaminases convert cytosine (C) to uracil (U) on DNA that is not shielded by nucleosomes or bound proteins. Accessible chromatin is therefore edited at a high rate, while protein-bound DNA is protected, leaving **footprints**. After PCR (which reads U as T) and sequencing, a deaminated cytosine appears as a C→T substitution relative to the reference; on the opposite strand the same event appears as G→A.

The signal is a mismatch: the reference carries C (or G), but the aligned read carries T (or A). The density of these mismatches reports accessibility, and local depletions reveal bound factors.

## Strand convention

The two mismatch patterns say **which strand of the molecule was deaminated** — not which way the observing read happened to align.

**Top-strand event.** A C on the plus strand is deaminated:

```
Reference (+):  5'—...C...—3'
                       ↓ deamination
Molecule (+):   5'—...T...—3'
```

Pattern: `ref_base == 'C'` and `read_base == 'T'`.

**Bottom-strand event.** A C on the minus strand is deaminated. That C pairs with a reference G, so in reference coordinates the event shows up at a G:

```
Reference (+):  5'—...G...—3'
Template (−):   3'—...C...—5'  ← deaminated here
```

Pattern: `ref_base == 'G'` and `read_base == 'A'`.

**Read orientation does not select between them.** Two facts make `is_reverse` the wrong test:

1. A BAM stores `SEQ` in **reference orientation** — a reverse-strand record has already been reverse-complemented — so both patterns are expressed against the plus strand regardless of how the read aligned.
2. Library prep fixes the deaminated U into a real **T:A** base pair, so after amplification *both* strands of the molecule carry the substitution and any read covering the position reports it.

A double-stranded deaminase edits both strands, so one fragment routinely carries both patterns at different positions. Every command therefore calls both patterns on every read. Measured on `data/ACCESS-ATAC/chr10.bam`, `is_reverse` predicts a read's dominant edit direction only about **65%** of the time; filtering on it discarded **35.5%** of editing events — events carrying the same `TC` enzyme-motif fingerprint as the ones kept, so genuine deamination rather than noise. `bam2fragment` used such a filter until 2026-09-10.

## Deamination-aware alignment (`index` + `align`)

Heavily edited reads align poorly to an unmodified reference, so DeamTools uses a bwa-meth-style three-letter strategy.

- **`index`** writes a *doubly-converted* reference: every chromosome appears twice, once with all C→T (prefixed `f`) and once with all G→A (prefixed `r`), then runs `bwa index` on it.
- **`align`** maps each read in **both** conversion directions and keeps whichever scores higher. Because the deaminase edits both strands, a single read can carry both `C→T` and `G→A` edits, so each read is emitted twice: single-end as `C→T` (`YC:Z:ct`) and `G→A` (`YC:Z:ga`); paired-end as two fragment orientations `f` = (R1 `C→T`, R2 `G→A`) and `r` = (R1 `G→A`, R2 `C→T`). The original sequence is stashed in `YS:Z:` and the candidate in `YC:Z:`, both carried through `bwa mem -C`.

After mapping, records are grouped by read name and the candidate with the higher primary alignment score (sum of the mates' `AS` for pairs) is kept; the original SEQ is restored from `YS`, the `f`/`r` prefix is stripped from RNAME/RNEXT, the `YS`/`YC` tags are dropped, and the BAM is sorted and indexed. Choosing one orientation per fragment keeps the mates on the same converted contig, so proper pairing is preserved.

### Read space vs reference space

`f` = (R1 `C→T`, R2 `G→A`) reads as though the two mates were converted
differently. **They are not.** Each candidate is *one* conversion applied to the
whole fragment in **reference space**; the labels differ only because R2's raw
FASTQ sequence is the reverse complement of the reference-space sequence, and
the conversion has to be written in the space the FASTQ is actually in.

A fragment whose **top** strand was deaminated, with R1 reading its left end and
R2 its right end:

```
reference (+)      ACGTCGATCGTTAGCCATGC
deaminated at      *      *      *
molecule, top      ATGTCGATTGTTAGCTATGC

R1 (read space)    ATGTCGATTG               ← already reference space
R2 (read space)              GCATAGCTAA     ← reverse complement of it
R2 in ref space              TTAGCTATGC
```

Candidate `f` means "convert this fragment `C→T` in reference space":

```
R1   read-space C→T   ATGTTGATTG
R2   read-space G→A   ACATAACTAA  ──bwa reverse-complements──▶  TTAGTTATGT
R2   ref-space  C→T                                             TTAGTTATGT   ✓ identical
```

The two agree because complementing turns `G→A` into `C→T`. Formally, for any
sequence `s`:

```
revcomp(s.replace("G", "A"))  ==  revcomp(s).replace("C", "T")
```

So `f` puts **both** mates on the `f<chrom>` contig and `r` puts both on
`r<chrom>`. That is the whole point of the label flip: it is what keeps the pair
on one contig.

### Why there are two candidates and not four

Converting each mate independently would give four combinations. The two extra
ones are not alternative hypotheses about the molecule — they place the mates on
*different* contigs, because a `C→T` mate goes to `f<chrom>` and a `G→A` mate to
`r<chrom>`. Measured on 300 simulated top-strand-deaminated fragments against
`data/test/chr20sim.deamtools.c2t`:

| candidate | mates on one contig | proper pairs | mean `AS(R1)+AS(R2)` |
|---|---|---|---|
| `f` = (R1 `ct`, R2 `ga`) | 278/279 | 273/279 | **200.0** |
| `r` = (R1 `ga`, R2 `ct`) | 272/279 | 226/279 | 159.2 |
| (R1 `ct`, R2 `ct`) | 20/279 | **0/279** | 177.4 |
| (R1 `ga`, R2 `ga`) | 19/279 | **0/279** | 176.9 |

Each mate still aligns somewhere on its own, so the invalid combinations collect
a respectable `AS` — 177.4 here, **above** the valid `r` candidate's 159.2. Adding
them to the take-best would therefore sometimes select a cross-contig,
non-proper pair over a correct one, and after the `f`/`r` prefix is stripped the
two mates would land on the same chromosome name carrying mate coordinates and
`TLEN` that mean nothing.

The reason only one choice exists per fragment is biological, not a
simplification: R1 and R2 report the *same molecule* in the same reference
orientation, so whatever edits it carries appear in both. The conversion is a
property of the fragment, so there is one choice to make, and take-best makes it.

What take-best cannot do is absorb **both** directions at once — a
double-stranded deaminase edits both strands, so a real fragment carries `C→T`
and `G→A` together and neither three-letter alphabet covers both. That is a
genuine limit of the three-letter approach, and it is why `deamtools` accuracy
declines gently with editing rate (99.26% → 99.11% from 0% to 100% editing)
rather than staying flat. Fixing *that* would need a different alignment model,
not more conversion candidates.

## Edit detection conventions

All commands iterate aligned positions with `pysam`'s `get_aligned_pairs(matches_only=True)`, so only matched (M/=/X) bases are compared — insertions, deletions, and clips are skipped, and query/reference bases are always directly comparable.

Two strand conventions are used, depending on the command:

| Command | Convention | Counted as an edit |
|---|---|---|
| `bam2bw`, `qc` | **strand-agnostic** | any reference `C→T` **or** `G→A` mismatch, regardless of read orientation |
| `bam2fragment` | **strand-agnostic, strand-resolved** | the same mismatches, reported in two separate columns so the edited strand is preserved |

Reads flagged unmapped, duplicate, QC-fail, secondary, or supplementary are always excluded, then `min_mapq` is applied per read; `min_baseq` gates individual bases.

### Mates are merged before counting

`bam2bw` and `qc` count per **fragment**, not per alignment record. Whenever the insert is shorter than twice the read length the two mates overlap, and every reference position in that overlap is reported by both of them — one position on one molecule, seen twice. Records are therefore grouped by read name and each fragment's mates are collapsed into one `reference position → base` map before anything is tallied, so an overlapping position contributes once.

Where the mates disagree at an overlap, the **higher base quality wins**: library prep turns the deaminated C into a real T:A pair that both strands carry, so the two mates are expected to agree and a disagreement means one of them misread.

This is not a rounding detail. On `data/ACCESS-ATAC/chr10.bam` the overlap accounted for **24.5%** of `qc`'s editing opportunities and **21.2%** of `bam2bw`'s count-mode signal. Because the overlap is the middle of the fragment rather than a random subset of positions, counting it twice biased every rate derived from it rather than merely inflating the totals.

Records whose mate never arrives — unpaired reads, a mate that was unmapped, filtered out, or placed on another contig — are counted as one-record fragments, so nothing is dropped. Single-end data is unaffected: one record is one fragment.

`bam2fragment` has always worked per fragment, and now shares the same merging helper.

## Signal generation (`bam2bw`)

For each region (whole chromosome or a merged BED interval) the reference is fetched once and reads are streamed via the BAM index.

- **count mode** — a per-base count of edits. With `--extend_size E > 0`, each event is broadcast symmetrically into a window of width `2E + 1` (clipped to the region), so the value at a base is the number of events within `E` bp.
- **ratio mode** — `edit_count / total_ACGT_coverage` at each base; positions whose total coverage is below `--min_coverage` are written as `0`. `--extend_size` is ignored in this mode.

Both are per fragment. The denominator is counted in the same pass as the numerator, over the same merged fragments, so the two cannot disagree about which reads they saw — and it now honours `--min_mapq`. Before 2026-09-10 the denominator came from pysam's `count_coverage`, which applies no mapping-quality filter, so ratio mode was dividing edits from MAPQ-passing reads by coverage that included reads the numerator had excluded.

BED intervals are merged before counting so a read spanning an overlap is not double-counted:

```
Input:   [0, 500)  [300, 800)  [1000, 1500)
Merged:  [0, 800)              [1000, 1500)
```

## Enzyme sequence bias

Deaminases have intrinsic flanking-sequence preferences (e.g. DddA strongly prefers `TC`). `qc` reports the per-trinucleotide edit fraction and renders a deaminase sequence-motif logo (the reference window around each editing event, centre excluded, `G→A` events reverse-complemented), both built directly from the BAM. A strongly skewed profile indicates that footprinting/occupancy analyses should correct for enzyme bias.

## Parallelism and output

Each chromosome (or region) is processed in its own thread via `concurrent.futures.ThreadPoolExecutor`; BAM and FASTA handles are opened independently per worker (pysam handles are not thread-safe). For BigWig output, regions are sorted into `(chrom, start, end)` order — as required by `pyBigWig` — and only non-zero bases are written (sparse, `span=1`); the header lists every chromosome even if it has no signal.

## Complexity

For *N* aligned reads of average length *L*, edit detection is O(*N* × *L*) per thread. `bam2bw` holds one float32 array per in-flight chromosome (≈ one float per base), so peak memory scales with the largest chromosomes times the number of threads.
