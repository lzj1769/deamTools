"""Sequence-to-edit (seq2edit) CNN modelling for DeamTools.

Learns a *DNA sequence -> per-base editing* map (the deaminase sequence bias)
and serves it for downstream bias correction.

Two architectures live here:

* :class:`~deamtools.seq2edit.bpnet.BPNet` in :mod:`deamtools.seq2edit.bpnet` --
  the ChromBPNet architecture (dilated residual stack, separate profile and
  counts heads). This is the model to use for new work.
* :class:`~deamtools.seq2edit.model.EditNet` in
  :mod:`deamtools.seq2edit.model` -- the earlier, much smaller net ported from
  ACCESS-ATAC's ``cnn_bias_model``, kept as a baseline. The sequence and signal
  encoding helpers also live in that module.

The workflow has three stages:

1. **train** -- fit a CNN on one-hot DNA windows against an editing-signal
   BigWig with a Poisson loss (implemented here; see :func:`run_train`).
2. **predict** -- score new sequences to produce an *expected* track *(planned)*.
3. **interpret** -- attribute predictions back to sequence *(planned)*.

The model and training loop require the optional ``torch`` dependency; install
the extra with ``pip install 'deamtools[seq2edit]'``.
"""

from deamtools.seq2edit.train import run_train

__all__ = ["run_train"]
