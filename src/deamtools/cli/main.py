from __future__ import annotations

import argparse

from deamtools.utils import get_version

from deamtools.preprocessing.bam2bw import run_bam2bw

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

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        help="Available subcommands (see 'deamtools <command> --help')",
    )

    _add_bam2bw_parser(subparsers)

    return parser


def _add_bam2bw_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "bam2bw",
        help=(
            "Convert a coordinate-sorted BAM file to a BigWig track of "
            "per-base C-to-T deamination counts."
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
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa --output sample.bw\n"
            "\n"
            "  # Region-restricted run with stricter quality filters and 4 threads\n"
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa \\\n"
            "      --regions peaks.bed --min_mapq 30 --min_baseq 30 \\\n"
            "      --threads 4 --output sample_peaks.bw\n"
            "\n"
            "  # Extend each editing site by 50 bp in both directions\n"
            "  deamtools bam2bw --bam sample.bam --fasta hg38.fa \\\n"
            "      --extend_size 50 --output sample_extended.bw\n"
            "\n"
            "notes:\n"
            "  * The BAM file must be coordinate-sorted and indexed (.bai).\n"
            "  * The FASTA file must be indexed with 'samtools faidx' (.fai).\n"
            "  * Output parent directories are created automatically if they do not exist."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        required=True,
        metavar="FILE",
        help=(
            "Path to the input BAM file. Must be coordinate-sorted and "
            "accompanied by an index file (.bai)."
        ),
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help=(
            "Path to the reference genome FASTA file used during alignment. "
            "Must be indexed with 'samtools faidx' (.fai)."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help=(
            "Path for the output BigWig file (.bw). "
            "Parent directories will be created if they do not exist."
        ),
    )
    parser.add_argument(
        "--chrom_sizes",
        required=True,
        metavar="FILE",
        help=(
            "Path to a chromosome sizes file (tab-delimited: chrom\\tsize). "
            "Chromosome sizes are required for BigWig creation."
        ),
    )
    parser.add_argument(
        "--regions",
        metavar="FILE",
        help=(
            "Path to a BED file defining genomic regions of interest. "
            "When provided, only reads overlapping these regions are processed, "
            "which can substantially reduce runtime for targeted analyses. "
            "If omitted, the entire genome is analysed."
        ),
    )
    parser.add_argument(
        "--extend_size",
        type=int,
        default=0,
        metavar="INT",
        help=(
            "Symmetrically extend each detected editing site by INT base pairs "
            "in both directions before writing to the BigWig. "
            "Useful for smoothing sparse signals or calling accessible windows "
            "around editing events. Default: %(default)s (no extension)."
        ),
    )
    parser.set_defaults(func=_run_bam2bw)


def _run_bam2bw(args: argparse.Namespace) -> int:
    run_bam2bw(
        bam_path=args.bam,
        fasta_path=args.fasta,
        bed_path=args.regions,
        output_path=args.output,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        threads=args.threads,
        log_level=args.log_level,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        return 1

    return int(args.func(args))