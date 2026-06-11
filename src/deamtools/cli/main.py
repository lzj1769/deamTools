from __future__ import annotations

import argparse
import logging
import shlex
import sys

from deamtools.align.align import run_align
from deamtools.align.index import run_index
from deamtools.motif.match import run_motif_matching
from deamtools.preprocessing.bam2bw import run_bam2bw
from deamtools.preprocessing.bam2fragment import run_bam2fragment
from deamtools.qc import run_qc
from deamtools.utils import get_version

logger = logging.getLogger(__name__)

# Argparse-internal attributes we don't want to print as user parameters.
_INTERNAL_ARG_KEYS = frozenset({"func", "command"})


def _log_invocation(args: argparse.Namespace) -> None:
    """Echo the invocation and resolved argument values, MACS2-style."""
    logger.info("# Command line: %s", " ".join(shlex.quote(a) for a in sys.argv))
    logger.info("# ARGUMENTS LIST:")
    for key, value in vars(args).items():
        if key in _INTERNAL_ARG_KEYS:
            continue
        logger.info("# %s = %s", key, value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="deamtools",
        description=(
            "DeamTools: a command-line toolkit for deamination-based chromatin\n"
            "accessibility analysis.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        metavar="LEVEL",
        help="Logging verbosity (DEBUG/INFO/WARNING/ERROR). Default: INFO.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        help="Available subcommands (see 'deamtools <command> --help')",
    )

    _add_index_parser(subparsers)
    _add_align_parser(subparsers)
    _add_bam2bw_parser(subparsers)
    _add_bam2fragment_parser(subparsers)
    _add_qc_parser(subparsers)
    _add_match_parser(subparsers)

    return parser


def _add_index_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "index",
        help="Build a deamination-aware BWA index for a reference FASTA.",
        description=(
            "Build a doubly-converted BWA index for deamination-aware alignment.\n"
            "\n"
            "The index contains two copies of every chromosome: one C-to-T converted\n"
            "(prefixed 'f') and one G-to-A converted (prefixed 'r'). This lets BWA-MEM\n"
            "map both top-strand-derived (C->T pattern) and bottom-strand-derived\n"
            "(G->A pattern in read orientation) deaminated reads against a single index.\n"
            "\n"
            "Outputs:\n"
            "  <fasta>.fai                            (always next to the FASTA;\n"
            "                                          required by the other commands)\n"
            "  <out_dir>/<out_name>.deamtools.c2t     (converted reference)\n"
            "  <out_dir>/<out_name>.deamtools.c2t.*   (BWA-MEM index files)\n"
            "\n"
            "--out_dir/--out_name default to the FASTA's directory and file name,\n"
            "so by default the index is written right next to the FASTA (the\n"
            "location 'deamtools align' looks in)."
        ),
        epilog=(
            "examples:\n"
            "  deamtools index --fasta hg38.fa\n"
            "  deamtools index --fasta hg38.fa --out_dir idx --out_name hg38\n"
            "  deamtools index --fasta hg38.fa --force\n"
            "\n"
            "notes:\n"
            "  * Requires 'bwa' and 'samtools' on PATH.\n"
            "  * Existing outputs are kept unless --force is given.\n"
            "  * The .fai is always written next to the FASTA."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Reference FASTA file to index.",
    )
    parser.add_argument(
        "--out_dir",
        metavar="DIR",
        help=(
            "Directory for the converted reference + BWA index. "
            "Default: the FASTA's directory."
        ),
    )
    parser.add_argument(
        "--out_name",
        metavar="NAME",
        help=(
            "Base name for the converted reference + BWA index. "
            "Default: the FASTA file name."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild even if existing index files are present.",
    )
    parser.set_defaults(func=_run_index)


def _add_align_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "align",
        help=(
            "Align deaminated reads (single- or paired-end) to a reference indexed "
            "with 'deamtools index'."
        ),
        description=(
            "Align deaminated sequencing reads to a deamtools-indexed reference.\n"
            "\n"
            "Reads are converted on the fly (read 1: C->T, read 2: G->A for paired-end)\n"
            "and aligned with BWA-MEM to the doubly-converted reference produced by\n"
            "'deamtools index'. Original read sequences are restored in the output BAM,\n"
            "which is sorted and indexed."
        ),
        epilog=(
            "examples:\n"
            "  # Paired-end\n"
            "  deamtools align --fasta hg38.fa --read1 r1.fq.gz --read2 r2.fq.gz \\\n"
            "      --out_dir results --out_name sample --threads 8\n"
            "\n"
            "  # Single-end with a read group\n"
            "  deamtools align --fasta hg38.fa --read1 reads.fq.gz \\\n"
            "      --read_group '@RG\\tID:s1\\tSM:sample1\\tLB:lib1\\tPL:ILLUMINA' \\\n"
            "      --out_dir results --out_name sample\n"
            "\n"
            "notes:\n"
            "  * Run 'deamtools index --fasta <ref>' once before aligning.\n"
            "  * Requires 'bwa' and 'samtools' on PATH."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Reference FASTA. Must have been indexed with 'deamtools index'.",
    )
    parser.add_argument(
        "--index",
        metavar="FILE",
        help=(
            "Path to the converted reference built by 'deamtools index' "
            "(<out_dir>/<out_name>.deamtools.c2t). Use this when the index was "
            "built with a custom --out_dir/--out_name. Default: next to the FASTA."
        ),
    )
    parser.add_argument(
        "--read1",
        required=True,
        metavar="FILE",
        help="FASTQ for read 1 (or for single-end reads). Plain or gzipped.",
    )
    parser.add_argument(
        "--read2",
        metavar="FILE",
        help="FASTQ for read 2 (paired-end). Omit for single-end alignment.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--out_name",
        required=True,
        metavar="NAME",
        help=(
            "Base name (without extension) for the output; writes a sorted, "
            "indexed <out_dir>/<out_name>.bam."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="INT",
        help="Total threads, split between bwa mem and samtools sort. Default: %(default)s.",
    )
    parser.add_argument(
        "--read_group",
        metavar="STR",
        help=(
            "Read group line passed to 'bwa mem -R', e.g. "
            "'@RG\\tID:s1\\tSM:sample1\\tLB:lib1\\tPL:ILLUMINA'."
        ),
    )
    parser.set_defaults(func=_run_align)


def _add_bam2bw_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "bam2bw",
        help=(
            "Convert a coordinate-sorted BAM file to a BigWig track of "
            "per-base C-to-T deamination counts or ratios."
        ),
        description=(
            "Convert aligned reads in BAM format to BigWig format by quantifying\n"
            "cytosine deamination (C-to-T editing) events at single-base resolution.\n"
            "\n"
            "For each read, the tool scans reference cytosine positions covered by the\n"
            "read and records positions where a C-to-U/T conversion is observed.\n"
            "The resulting per-base editing counts are written to a BigWig file suitable\n"
            "for genome-browser visualisation or downstream signal analysis."
        ),
        epilog=(
            "examples:\n"
            "  # Whole-genome run with default quality thresholds\n"
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa \\\n"
            "      --out_dir results --out_name sample\n"
            "\n"
            "  # Region-restricted run with stricter quality filters and 4 threads\n"
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa \\\n"
            "      --regions peaks.bed --min_mapq 30 --min_baseq 30 \\\n"
            "      --threads 4 --out_dir results --out_name sample_peaks\n"
            "\n"
            "  # Extend each editing site by 50 bp in both directions\n"
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa \\\n"
            "      --extend_size 50 --out_dir results --out_name sample_extended\n"
            "\n"
            "notes:\n"
            "  * The BAM file must be coordinate-sorted and indexed (.bai).\n"
            "  * The FASTA file must be indexed with 'samtools faidx' (.fai).\n"
            "  * Output parent directories are created automatically if they do not exist."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required inputs
    parser.add_argument(
        "--bam",
        required=True,
        metavar="FILE",
        help="Path to the coordinate-sorted, indexed BAM file (.bai required).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Path to the reference FASTA file indexed with 'samtools faidx' (.fai required).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--out_name",
        required=True,
        metavar="NAME",
        help=(
            "Base name (without extension) for the output BigWig; writes "
            "<out_dir>/<out_name>.bw."
        ),
    )

    # Optional inputs
    parser.add_argument(
        "--chrom_sizes",
        metavar="FILE",
        help=(
            "Tab-delimited chromosome sizes file (chrom\\tsize). "
            "If omitted, sizes are inferred from the BAM header."
        ),
    )
    parser.add_argument(
        "--regions",
        metavar="FILE",
        help=(
            "BED file of genomic regions to restrict analysis to. "
            "When omitted, the entire genome is processed."
        ),
    )

    # Signal options
    parser.add_argument(
        "--mode",
        default="count",
        choices=["count", "ratio"],
        metavar="MODE",
        help=(
            "Signal to write to BigWig. 'count' (default) writes the raw number "
            "of deamination events at each base (any C->T or G->A reference "
            "mismatch). 'ratio' writes edit_count / total_coverage at each base, "
            "where total_coverage is the sum of ACGT read coverage and positions "
            "below --min_coverage are reported as 0."
        ),
    )
    parser.add_argument(
        "--extend_size",
        type=int,
        default=0,
        metavar="INT",
        help=(
            "Symmetrically extend each detected editing site by INT base pairs "
            "before writing to BigWig. Only applied in --mode count; ignored in "
            "--mode ratio. Default: %(default)s (no extension)."
        ),
    )
    parser.add_argument(
        "--min_coverage",
        type=int,
        default=10,
        metavar="INT",
        help=(
            "Minimum total ACGT coverage required to report a non-zero ratio in "
            "--mode ratio. Positions below this threshold are written as 0. "
            "Ignored in --mode count. Default: %(default)s."
        ),
    )

    # Quality filters
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        metavar="INT",
        help="Minimum read mapping quality to include. Default: %(default)s.",
    )
    parser.add_argument(
        "--min_baseq",
        type=int,
        default=20,
        metavar="INT",
        help="Minimum base quality at a position to count an editing event. Default: %(default)s.",
    )

    # Performance
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads for parallel region processing. Default: %(default)s.",
    )

    parser.set_defaults(func=_run_bam2bw)


def _add_bam2fragment_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "bam2fragment",
        help=(
            "Convert a coordinate-sorted BAM file to a per-fragment editing-signal "
            "table."
        ),
        description=(
            "Convert aligned reads in BAM format to a tab-delimited fragment table.\n"
            "\n"
            "Each row reports a unique fragment defined by (chrom, start, end,\n"
            "editing positions [, barcode]). The 'count' column gives the number\n"
            "of reads/pairs producing that exact signature; the 'edits' column is\n"
            "a '|'-separated list of 0-based reference positions where a C->T\n"
            "(forward read) or G->A (reverse read) deamination event was observed.\n"
            "\n"
            "Fragments are formed by pairing properly-paired reads (start = min of\n"
            "the two read starts, end = max of the two read ends); for unpaired\n"
            "BAMs each read is treated as a single-end fragment."
        ),
        epilog=(
            "examples:\n"
            "  # Bulk fragment table (no barcode)\n"
            "  deamtools bam2fragment --bam sample.bam --fasta hg38.fa \\\n"
            "      --out_dir results --out_name sample\n"
            "\n"
            "  # Single-cell, gzip-compressed, with 10x-style barcode tag\n"
            "  deamtools bam2fragment --bam sample.bam --fasta hg38.fa \\\n"
            "      --barcode --barcode_tag CB --gzip \\\n"
            "      --out_dir results --out_name sample\n"
            "\n"
            "notes:\n"
            "  * The BAM must be coordinate-sorted and indexed (.bai).\n"
            "  * The FASTA must be indexed with 'samtools faidx' (.fai).\n"
            "  * Writes <out_dir>/<out_name>.tsv (or .tsv.gz with --gzip)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        required=True,
        metavar="FILE",
        help="Coordinate-sorted, indexed BAM file (.bai required).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Reference FASTA file indexed with 'samtools faidx' (.fai required).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--out_name",
        required=True,
        metavar="NAME",
        help=(
            "Base name (without extension) for the output; writes "
            "<out_dir>/<out_name>.tsv (or .tsv.gz with --gzip)."
        ),
    )
    parser.add_argument(
        "--gzip",
        action="store_true",
        help="Write the fragment table gzip-compressed (<out_name>.tsv.gz).",
    )

    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        metavar="INT",
        help="Minimum read mapping quality. Default: %(default)s.",
    )
    parser.add_argument(
        "--min_baseq",
        type=int,
        default=20,
        metavar="INT",
        help="Minimum base quality at a position to count an editing event. Default: %(default)s.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads for parallel chromosome processing. Default: %(default)s.",
    )

    parser.add_argument(
        "--barcode",
        action="store_true",
        help=(
            "Include a barcode column in the output. With this flag the output "
            "format becomes 'chrom\\tstart\\tend\\tbarcode\\tcount\\tedits' "
            "(10x fragments-style). Fragments without the barcode tag are written "
            "with '.' as the barcode."
        ),
    )
    parser.add_argument(
        "--barcode_tag",
        default="CB",
        metavar="TAG",
        help="BAM tag carrying the cell barcode. Default: %(default)s (10x convention).",
    )

    parser.set_defaults(func=_run_bam2fragment)


def _add_qc_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "qc",
        help=(
            "Compute quality-control metrics (editing rate, enzyme context "
            "bias, fragment sizes, TSS enrichment) for a deaminase BAM."
        ),
        description=(
            "Summarise a coordinate-sorted BAM and its reference FASTA into the\n"
            "quality-control metrics most useful for a deaminase footprinting\n"
            "experiment:\n"
            "\n"
            "  * Read statistics: totals plus duplicate, properly-paired,\n"
            "    secondary and supplementary fractions.\n"
            "  * Editing statistics: genome-wide deamination rate (edits over\n"
            "    editable C/G opportunities) and edits-per-read distribution.\n"
            "  * Trinucleotide context bias: edit fraction per cytosine-centred\n"
            "    trinucleotide (G->A events reverse-complemented to the C->T\n"
            "    orientation) -- the enzyme's sequence-preference fingerprint.\n"
            "  * Fragment-length distribution from properly-paired reads.\n"
            "  * TSS enrichment (optional, when --tss is supplied).\n"
            "\n"
            "Two files are written: a machine-readable JSON and a self-contained,\n"
            "MultiQC-style HTML report that embeds the summary figure and\n"
            "documents the meaning of every metric inline\n"
            "(<out_dir>/<out_name>.json and .html)."
        ),
        epilog=(
            "examples:\n"
            "  # Core metrics from BAM + FASTA\n"
            "  deamtools qc --bam sample.bam --fasta hg38.fa \\\n"
            "      --out_dir results --out_name sample\n"
            "\n"
            "  # Add TSS enrichment and run on 4 threads\n"
            "  deamtools qc --bam sample.bam --fasta hg38.fa --tss tss.bed \\\n"
            "      --threads 4 --out_dir results --out_name sample\n"
            "\n"
            "notes:\n"
            "  * The BAM must be coordinate-sorted and indexed (.bai).\n"
            "  * The FASTA must be indexed with 'samtools faidx' (.fai).\n"
            "  * The TSS BED is read as (chrom, start, end); the TSS is taken as\n"
            "    the interval midpoint."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        required=True,
        metavar="FILE",
        help="Coordinate-sorted, indexed BAM file (.bai required).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Reference FASTA file indexed with 'samtools faidx' (.fai required).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--out_name",
        required=True,
        metavar="NAME",
        help=(
            "Base name (without extension) for the outputs; writes "
            "<out_dir>/<out_name>.json and <out_dir>/<out_name>.html."
        ),
    )
    parser.add_argument(
        "--tss",
        metavar="FILE",
        help=(
            "BED file of transcription start sites. When supplied, an "
            "ATAC-style TSS enrichment score and profile are computed."
        ),
    )
    parser.add_argument(
        "--tss_flank",
        type=int,
        default=2000,
        metavar="INT",
        help="Half-width in bp of the window around each TSS. Default: %(default)s.",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        metavar="INT",
        help="Minimum read mapping quality. Default: %(default)s.",
    )
    parser.add_argument(
        "--min_baseq",
        type=int,
        default=20,
        metavar="INT",
        help=(
            "Minimum base quality for a position to count as an editing "
            "opportunity. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads for parallel chromosome processing. Default: %(default)s.",
    )
    parser.add_argument(
        "--no_plot",
        action="store_true",
        help="Skip rendering/embedding the summary figure in the HTML report.",
    )

    parser.set_defaults(func=_run_qc)


def _add_match_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "match",
        help=(
            "Scan genomic regions for transcription-factor motif matches and "
            "write a BED of binding sites."
        ),
        description=(
            "Scan the reference sequence of a set of regions (e.g. peaks) for\n"
            "transcription-factor motif occurrences using MOODS, and write the\n"
            "hits as a BED of motif-predicted binding sites.\n"
            "\n"
            "Each motif's count matrix is converted to a log-odds matrix against\n"
            "a flat background; a per-motif score threshold is derived from\n"
            "--p_value and both strands are scanned. Motifs are fetched from\n"
            "JASPAR (requires the optional 'pyjaspar' package)."
        ),
        epilog=(
            "examples:\n"
            "  deamtools match --fasta hg38.fa --regions peaks.bed \\\n"
            "      --out_dir results --out_name mpbs\n"
            "\n"
            "  deamtools match --fasta hg38.fa --regions peaks.bed \\\n"
            "      --collection CORE --tax_group vertebrates --p_value 1e-4 \\\n"
            "      --out_dir results --out_name mpbs\n"
            "\n"
            "notes:\n"
            "  * The FASTA must be indexed with 'samtools faidx' (.fai).\n"
            "  * Motif fetching needs 'pyjaspar' (pip install pyjaspar).\n"
            "  * Output is 6-column BED: chrom, start, end, motif, score, strand."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Reference FASTA indexed with 'samtools faidx' (.fai required).",
    )
    parser.add_argument(
        "--regions",
        required=True,
        metavar="FILE",
        help="BED file of regions to scan (overlapping intervals are merged).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "--out_name",
        required=True,
        metavar="NAME",
        help=(
            "Base name (without extension) for the output; writes "
            "<out_dir>/<out_name>.bed."
        ),
    )
    parser.add_argument(
        "--jaspar_release",
        default="JASPAR2024",
        metavar="STR",
        help="JASPAR release to fetch motifs from. Default: %(default)s.",
    )
    parser.add_argument(
        "--collection",
        default="CORE",
        metavar="STR",
        help="JASPAR motif collection. Default: %(default)s.",
    )
    parser.add_argument(
        "--tax_group",
        nargs="*",
        metavar="GROUP",
        help="JASPAR taxonomic group(s). Default: vertebrates.",
    )
    parser.add_argument(
        "--p_value",
        type=float,
        default=1e-4,
        metavar="FLOAT",
        help="Significance threshold for motif hits. Default: %(default)s.",
    )
    parser.set_defaults(func=_run_match)


def _run_index(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_index(
        fasta_path=args.fasta,
        out_dir=args.out_dir,
        out_name=args.out_name,
        force=args.force,
    )
    return 0


def _run_align(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_align(
        fasta_path=args.fasta,
        read1=args.read1,
        read2=args.read2,
        out_dir=args.out_dir,
        out_name=args.out_name,
        threads=args.threads,
        read_group=args.read_group,
        index_path=args.index,
    )
    return 0


def _run_bam2bw(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_bam2bw(
        bam_path=args.bam,
        fasta_path=args.fasta,
        out_dir=args.out_dir,
        out_name=args.out_name,
        chrom_sizes_path=args.chrom_sizes,
        bed_path=args.regions,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        extend_size=args.extend_size,
        threads=args.threads,
        mode=args.mode,
        min_coverage=args.min_coverage,
    )
    return 0


def _run_bam2fragment(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_bam2fragment(
        bam_path=args.bam,
        fasta_path=args.fasta,
        out_dir=args.out_dir,
        out_name=args.out_name,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        threads=args.threads,
        barcode=args.barcode,
        barcode_tag=args.barcode_tag,
        gzip=args.gzip,
    )
    return 0


def _run_qc(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_qc(
        bam_path=args.bam,
        fasta_path=args.fasta,
        out_dir=args.out_dir,
        out_name=args.out_name,
        tss_path=args.tss,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        threads=args.threads,
        tss_flank=args.tss_flank,
        plot=not args.no_plot,
    )
    return 0


def _run_match(args: argparse.Namespace) -> int:
    _log_invocation(args)
    run_motif_matching(
        fasta_path=args.fasta,
        bed_path=args.regions,
        out_dir=args.out_dir,
        out_name=args.out_name,
        release=args.jaspar_release,
        collection=args.collection,
        tax_group=args.tax_group,
        p_value=args.p_value,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=getattr(logging, args.log_level),
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if not hasattr(args, "func"):
        parser.print_help()
        return 1

    return int(args.func(args))
