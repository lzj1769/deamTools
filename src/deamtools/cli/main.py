from __future__ import annotations

import argparse
import logging

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
        "--output",
        required=True,
        metavar="FILE",
        help="Output BigWig file path (.bw). Parent directories are created automatically.",
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
        "--extend_size",
        type=int,
        default=0,
        metavar="INT",
        help=(
            "Symmetrically extend each detected editing site by INT base pairs "
            "before writing to BigWig. Default: %(default)s (no extension)."
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
        help="Number of threads for parallel chromosome processing. Default: %(default)s.",
    )

    parser.set_defaults(func=_run_bam2bw)


def _run_bam2bw(args: argparse.Namespace) -> int:
    run_bam2bw(
        bam_path=args.bam,
        fasta_path=args.fasta,
        output_path=args.output,
        chrom_sizes_path=args.chrom_sizes,
        bed_path=args.regions,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        extend_size=args.extend_size,
        threads=args.threads,
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
