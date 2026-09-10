"""A ChromBPNet-architecture model for per-base editing prediction.

:class:`BPNet` is a PyTorch port of the architecture ChromBPNet uses
(``chrombpnet/training/models/bpnet_model.py`` in kundajelab/chrombpnet, itself
derived from BPNet, Avsec et al. 2021): one wide convolution, a stack of dilated
convolutions with residual connections, and **two output heads**.

The two heads factor the prediction into *where* and *how much*:

* the **profile head** predicts logits over positions, which softmax to the
  *shape* of the signal within the window, trained against a multinomial
  likelihood; and
* the **counts head** predicts a single scalar, the log of the *total* signal in
  the window, trained against MSE.

That factorisation is why the architecture suits deaminase bias modelling. The
enzyme's sequence preference is a statement about which cytosines get edited
relative to their neighbours, which is exactly the profile; total editing in a
window is dominated by accessibility and copy number, which is exactly what the
counts head absorbs and keeps out of the profile.

Everything here is valid-padded, so the output window is strictly narrower than
the input: the model sees flanking context it never predicts on. Use
:func:`required_input_len` to size the input, or :func:`profile_out_len` to find
what a given input yields; :class:`BPNet` validates the combination and reports
the arithmetic when it does not work.

Requires the optional ``torch`` dependency: ``pip install 'deamtools[seq2edit]'``.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import Tensor, nn

__all__ = [
    "BPNet",
    "BPNetLoss",
    "body_out_len",
    "count_loss",
    "estimate_counts_loss_weight",
    "multinomial_nll",
    "profile_out_len",
    "required_input_len",
]


# --------------------------------------------------------------------------- #
# Length arithmetic.
#
# Every convolution is valid-padded, so each one shortens the sequence by a
# fixed amount that depends only on the hyper-parameters. Getting this wrong is
# the usual way to waste a training run, so it is computed up front and checked
# in the constructor rather than discovered at the first forward pass.
# --------------------------------------------------------------------------- #
def body_out_len(
    input_len: int, conv1_kernel_size: int = 21, n_dil_layers: int = 8
) -> int:
    """Length after the first convolution and the dilated stack.

    The first convolution removes ``conv1_kernel_size - 1`` positions. Dilated
    layer ``i`` (kernel 3, dilation ``2**i``) has an effective width of
    ``1 + 2 * 2**i`` and so removes ``2 * 2**i``; over ``i = 1..n`` that sums to
    ``2**(n + 2) - 4``.
    """
    return input_len - (conv1_kernel_size - 1) - (2 ** (n_dil_layers + 2) - 4)


def profile_out_len(
    input_len: int,
    conv1_kernel_size: int = 21,
    n_dil_layers: int = 8,
    profile_kernel_size: int = 75,
) -> int:
    """Widest output window an input of ``input_len`` can produce.

    This is the length before the profile head's final centre crop, i.e. the
    maximum usable ``output_len``.
    """
    return body_out_len(input_len, conv1_kernel_size, n_dil_layers) - (
        profile_kernel_size - 1
    )


def required_input_len(
    output_len: int,
    conv1_kernel_size: int = 21,
    n_dil_layers: int = 8,
    profile_kernel_size: int = 75,
) -> int:
    """Smallest input length that yields exactly ``output_len``, with no crop.

    Inverse of :func:`profile_out_len`. With ChromBPNet's defaults an
    ``output_len`` of 1000 requires an input of 2114.
    """
    trim = (
        (conv1_kernel_size - 1)
        + (2 ** (n_dil_layers + 2) - 4)
        + (profile_kernel_size - 1)
    )
    return output_len + trim


def _centre_crop(x: Tensor, target_len: int) -> Tensor:
    """Crop the last dimension of ``x`` to ``target_len``, keeping the centre."""
    excess = x.shape[-1] - target_len
    if excess == 0:
        return x
    if excess < 0:
        raise ValueError(f"cannot crop length {x.shape[-1]} up to {target_len}")
    left = excess // 2
    return x[..., left : left + target_len]


# --------------------------------------------------------------------------- #
# The network.
# --------------------------------------------------------------------------- #
class BPNet(nn.Module):
    """ChromBPNet-architecture network: dilated residual body, two heads.

    Consumes a one-hot window of shape ``(batch, input_len, 4)`` and returns a
    ``(profile_logits, log_counts)`` pair, shaped ``(batch, n_tasks,
    output_len)`` and ``(batch, n_tasks)``.

    The profile output is **logits, not rates**. Softmax over the position axis
    gives the predicted shape; multiply by ``expm1(log_counts)`` to recover an
    expected per-base signal (:meth:`predict_signal` does this).

    Parameters
    ----------
    input_len : int
        Input window width in bp. Must be wide enough to leave ``output_len``
        positions after the valid-padded stack; see :func:`required_input_len`.
    output_len : int
        Width of the predicted profile in bp, centred inside the input window.
    filters : int
        Channel width of the first convolution and of every dilated layer.
    n_dil_layers : int
        Number of dilated residual layers; dilation doubles from 2 up to
        ``2**n_dil_layers``.
    conv1_kernel_size : int
        Kernel width of the first convolution.
    profile_kernel_size : int
        Kernel width of the profile head's convolution.
    n_tasks : int
        Number of output tracks. One per editing signal being modelled.

    Raises
    ------
    ValueError
        If ``input_len`` is too short for ``output_len``, or if the leftover
        crop is odd (which would make the output window off-centre by half a
        base). The message reports the arithmetic and the nearest usable input
        length.
    """

    def __init__(
        self,
        input_len: int = 2114,
        output_len: int = 1000,
        filters: int = 512,
        n_dil_layers: int = 8,
        conv1_kernel_size: int = 21,
        profile_kernel_size: int = 75,
        n_tasks: int = 1,
    ) -> None:
        super().__init__()

        pre_crop = profile_out_len(
            input_len, conv1_kernel_size, n_dil_layers, profile_kernel_size
        )
        if pre_crop < output_len:
            needed = required_input_len(
                output_len, conv1_kernel_size, n_dil_layers, profile_kernel_size
            )
            raise ValueError(
                f"input_len={input_len} yields at most {pre_crop} output "
                f"positions, which is short of output_len={output_len}. "
                f"Use input_len={needed} (or larger by an even amount)."
            )
        if (pre_crop - output_len) % 2 != 0:
            raise ValueError(
                f"input_len={input_len} leaves {pre_crop} positions for an "
                f"output_len={output_len} window, an odd difference, so the "
                "window cannot be centred. Change input_len by one."
            )

        self.input_len = input_len
        self.output_len = output_len
        self.filters = filters
        self.n_dil_layers = n_dil_layers
        self.conv1_kernel_size = conv1_kernel_size
        self.profile_kernel_size = profile_kernel_size
        self.n_tasks = n_tasks

        self.conv1 = nn.Conv1d(4, filters, kernel_size=conv1_kernel_size)
        self.dilated = nn.ModuleList(
            nn.Conv1d(filters, filters, kernel_size=3, dilation=2**i)
            for i in range(1, n_dil_layers + 1)
        )
        self.profile_conv = nn.Conv1d(filters, n_tasks, kernel_size=profile_kernel_size)
        self.counts_dense = nn.Linear(filters, n_tasks)

    def forward(self, x: Tensor) -> tuple[Tensor, Tensor]:
        """Return ``(profile_logits, log_counts)`` for a one-hot batch.

        Parameters
        ----------
        x : torch.Tensor
            One-hot input, ``(batch, input_len, 4)``.

        Returns
        -------
        profile_logits : torch.Tensor
            ``(batch, n_tasks, output_len)``. Logits over positions, not rates.
        log_counts : torch.Tensor
            ``(batch, n_tasks)``. Predicted ``log1p`` of the window total.
        """
        # (batch, len, 4) -> (batch, 4, len) for Conv1d.
        x = x.permute(0, 2, 1)
        x = F.relu(self.conv1(x))

        # Each dilated layer is narrower than its input, so the residual branch
        # is centre-cropped to match before the add.
        for conv in self.dilated:
            conv_x = F.relu(conv(x))
            x = _centre_crop(x, conv_x.shape[-1]) + conv_x

        profile = _centre_crop(self.profile_conv(x), self.output_len)
        log_counts = self.counts_dense(x.mean(dim=-1))
        return profile, log_counts

    @torch.no_grad()
    def predict_signal(self, x: Tensor) -> Tensor:
        """Expected per-base signal, ``(batch, n_tasks, output_len)``.

        Combines the two heads the way the training objective factors them:
        the softmaxed profile is a distribution over positions, scaled by the
        window total recovered from the counts head. This is the track that
        downstream bias correction divides out.
        """
        profile, log_counts = self(x)
        shape = F.softmax(profile, dim=-1)
        total = torch.expm1(log_counts).clamp_min(0.0).unsqueeze(-1)
        return shape * total

    @property
    def config(self) -> dict:
        """Hyper-parameters needed to rebuild this model for inference."""
        return {
            "input_len": self.input_len,
            "output_len": self.output_len,
            "filters": self.filters,
            "n_dil_layers": self.n_dil_layers,
            "conv1_kernel_size": self.conv1_kernel_size,
            "profile_kernel_size": self.profile_kernel_size,
            "n_tasks": self.n_tasks,
        }


# --------------------------------------------------------------------------- #
# Losses.
# --------------------------------------------------------------------------- #
def multinomial_nll(logits: Tensor, true_counts: Tensor) -> Tensor:
    """Multinomial negative log-likelihood of ``true_counts`` under ``logits``.

    The profile head is scored on *shape only*: the observed window total is
    taken as the multinomial's ``n``, so the loss cannot be reduced by getting
    the overall level right. That is the counts head's job, which is what keeps
    the two heads from learning the same thing.

    Parameters
    ----------
    logits : torch.Tensor
        Predicted logits over positions, ``(..., positions)``.
    true_counts : torch.Tensor
        Observed per-base counts, same shape as ``logits``.

    Returns
    -------
    torch.Tensor
        Scalar mean NLL over all leading dimensions.
    """
    log_probs = F.log_softmax(logits, dim=-1)
    total = true_counts.sum(dim=-1)
    # Multinomial coefficient: constant in the parameters, kept so the reported
    # value is a true log-likelihood and is comparable with other runs.
    log_coeff = torch.lgamma(total + 1) - torch.lgamma(true_counts + 1).sum(dim=-1)
    log_lik = log_coeff + (true_counts * log_probs).sum(dim=-1)
    return -log_lik.mean()


def count_loss(log_counts: Tensor, true_counts: Tensor) -> Tensor:
    """MSE between predicted log-counts and ``log1p`` of the observed total.

    Parameters
    ----------
    log_counts : torch.Tensor
        Counts-head output, ``(batch, n_tasks)``.
    true_counts : torch.Tensor
        Observed per-base counts, ``(batch, n_tasks, positions)``.
    """
    return F.mse_loss(log_counts, torch.log1p(true_counts.sum(dim=-1)))


def estimate_counts_loss_weight(true_counts: Tensor, scale: float = 10.0) -> float:
    """A starting ``counts_weight`` for :class:`BPNetLoss`.

    The two losses are on unrelated scales: the multinomial NLL grows with the
    number of reads in a window while the counts MSE is order 1. ChromBPNet
    handles this by weighting the counts term by roughly the median window total
    over ``scale``; the same heuristic is reproduced here. It is a starting
    point, not a tuned value.

    Parameters
    ----------
    true_counts : torch.Tensor
        Observed per-base counts over a representative sample of windows,
        ``(n_windows, ..., positions)``.
    scale : float
        Divisor applied to the median window total.
    """
    totals = true_counts.sum(dim=-1).flatten().float()
    return float(totals.median().item() / scale)


class BPNetLoss(nn.Module):
    """Combined profile + counts objective.

    ``loss = multinomial_nll(profile) + counts_weight * mse(log_counts)``.

    ``counts_weight`` matters: leave it at 1 and the multinomial term, which
    scales with read depth, drowns the counts term out. Size it with
    :func:`estimate_counts_loss_weight`.

    Parameters
    ----------
    counts_weight : float
        Weight on the counts term.
    """

    def __init__(self, counts_weight: float = 1.0) -> None:
        super().__init__()
        self.counts_weight = counts_weight

    def forward(
        self,
        profile_logits: Tensor,
        log_counts: Tensor,
        true_counts: Tensor,
    ) -> tuple[Tensor, Tensor, Tensor]:
        """Return ``(total, profile_term, counts_term)``.

        The unweighted components are returned alongside the total so training
        can log them separately; the two move on very different scales and a
        single number hides which one is still improving.
        """
        profile_term = multinomial_nll(profile_logits, true_counts)
        counts_term = count_loss(log_counts, true_counts)
        total = profile_term + self.counts_weight * counts_term
        return total, profile_term, counts_term
