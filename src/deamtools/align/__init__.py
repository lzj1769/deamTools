"""Deamination-aware alignment using a bwa-meth-style strategy."""

from deamtools.align.align import run_align
from deamtools.align.index import run_index

__all__ = ["run_align", "run_index"]
