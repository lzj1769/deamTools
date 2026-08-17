"""Sequence-to-edit (seq2edit) CNN modelling for DeamTools.

Learns a *DNA sequence -> per-base editing* map (the deaminase sequence bias)
and serves it for downstream bias correction. The :class:`EditNet` model and its
encoding helpers live in :mod:`deamtools.seq2edit.model`. The workflow has three
stages:

1. **train** -- fit a CNN on one-hot DNA windows against an editing-signal
   BigWig with a Poisson loss (implemented here; see :func:`run_train`).
2. **predict** -- score new sequences to produce an *expected* track *(planned)*.
3. **interpret** -- attribute predictions back to sequence *(planned)*.

The model and training loop require the optional ``torch`` dependency; install
the extra with ``pip install 'deamtools[seq2edit]'``.
"""

from deamtools.seq2edit.train import run_train

__all__ = ["run_train"]
