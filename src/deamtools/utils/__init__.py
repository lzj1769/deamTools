"""Utility helpers for DeamTools."""

from deamtools.utils._logging import logger
from deamtools.utils.chromosome import (
    get_chrom_sizes_from_bam,
    get_chrom_sizes_from_file,
)
from deamtools.utils.fragments import (
    Fragment,
    iter_fragments,
    mates_can_pair,
    merge_fragment_bases,
)
from deamtools.utils.parallel import run_jobs
from deamtools.utils.regions import BED_COLUMNS, _load_regions
from deamtools.utils.version import get_version

__all__ = [
    "get_version",
    "logger",
    "get_chrom_sizes_from_bam",
    "get_chrom_sizes_from_file",
    "BED_COLUMNS",
    "_load_regions",
    "Fragment",
    "iter_fragments",
    "mates_can_pair",
    "merge_fragment_bases",
    "run_jobs",
]
