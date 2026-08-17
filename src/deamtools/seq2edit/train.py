"""Training loop for the seq2edit CNN model.

Trains :class:`deamtools.seq2edit.model.EditNet` to predict per-base editing
signal from one-hot DNA, following the ACCESS-ATAC ``cnn_bias_model`` recipe but
with a **Poisson** negative-log-likelihood loss instead of MSE: editing signal
is count-like, so the model outputs a positive per-base rate ``lambda`` and is
fit to maximise the Poisson likelihood of the observed counts. Optimisation uses
Adam (lr ``3e-4``, weight decay ``1e-4``) and a ``ReduceLROnPlateau`` schedule.
The best-validation checkpoint is written to ``<out_dir>/<out_name>.pth``
together with the loss history and the model config needed to rebuild it for
prediction/interpretation.

Requires the optional ``torch`` dependency: install with
``pip install 'deamtools[seq2edit]'``.
"""

from __future__ import annotations

import logging
import os

import numpy as np

logger = logging.getLogger(__name__)

_TORCH_HINT = (
    "PyTorch is required for 'deamtools seq2edit'. Install the optional extra "
    "with: pip install 'deamtools[seq2edit]'"
)


def _require_torch():
    """Import torch with an actionable error when the extra is not installed."""
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - exercised only without torch
        raise ImportError(_TORCH_HINT) from exc
    return torch


def _select_device(device: str | None):
    """Resolve the compute device, preferring CUDA, then Apple MPS, then CPU."""
    import torch

    if device is not None:
        return torch.device(device)
    if torch.cuda.is_available():
        return torch.device("cuda")
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _split_train_valid(
    x: np.ndarray, y: np.ndarray, valid_fraction: float, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Random hold-out split used when no explicit validation regions are given."""
    rng = np.random.default_rng(seed)
    n = x.shape[0]
    perm = rng.permutation(n)
    n_valid = max(1, int(round(n * valid_fraction)))
    valid_idx, train_idx = perm[:n_valid], perm[n_valid:]
    return x[train_idx], y[train_idx], x[valid_idx], y[valid_idx]


def _run_epoch(model, loader, loss_fn, device, optimizer=None) -> float:
    """One pass over ``loader``; trains when ``optimizer`` is given, else evals."""
    import torch

    train = optimizer is not None
    model.train(train)
    total, n = 0.0, 0
    with torch.set_grad_enabled(train):
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)
            loss = loss_fn(pred, yb)
            if train:
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
            total += loss.item() * xb.shape[0]
            n += xb.shape[0]
    return total / max(n, 1)


def run_train(
    fasta_path: str,
    bigwig_path: str,
    train_regions: str,
    out_dir: str,
    out_name: str,
    valid_regions: str | None = None,
    seq_len: int = 128,
    step: int | None = None,
    n_filters: int = 32,
    kernel_size: int = 5,
    epochs: int = 200,
    batch_size: int = 512,
    lr: float = 3e-4,
    weight_decay: float = 1e-4,
    lr_patience: int = 10,
    min_lr: float = 1e-5,
    valid_fraction: float = 0.1,
    num_workers: int = 0,
    seed: int = 42,
    device: str | None = None,
) -> str:
    """Train the seq2edit bias model and write the best checkpoint.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA indexed with ``samtools faidx``.
    bigwig_path : str
        Per-base editing-signal BigWig (ideally a naked/deproteinised-DNA
        control, so the model learns enzyme sequence bias).
    train_regions : str
        BED file of regions to draw training windows from.
    out_dir, out_name : str
        Checkpoint is written to ``<out_dir>/<out_name>.pth``.
    valid_regions : str, optional
        BED file of validation regions. When omitted, a ``valid_fraction``
        random hold-out of the training windows is used instead.
    seq_len, step, n_filters, kernel_size : int
        Data/architecture sizes (see :class:`EditNet` and
        :func:`build_dataset`).
    epochs, batch_size, lr, weight_decay, lr_patience, min_lr : numeric
        Optimisation hyper-parameters (ACCESS-ATAC defaults).
    valid_fraction : float
        Hold-out fraction when ``valid_regions`` is not supplied.
    num_workers, seed, device : misc
        DataLoader workers, RNG seed, and compute device override
        (``"cpu"``/``"cuda"``/``"mps"``; auto-detected when ``None``).

    Returns
    -------
    str
        Path to the written ``.pth`` checkpoint.
    """
    torch = _require_torch()
    from torch import nn

    from deamtools.seq2edit.dataset import DNASeqDataset, get_dataloader
    from deamtools.seq2edit.model import EditNet, build_dataset

    torch.manual_seed(seed)
    np.random.seed(seed)
    os.makedirs(out_dir, exist_ok=True)

    # --- data ---------------------------------------------------------------
    logger.info("Building training windows from %s", train_regions)
    x_train, y_train, _ = build_dataset(
        fasta_path, bigwig_path, train_regions, seq_len=seq_len, step=step
    )
    if valid_regions is not None:
        logger.info("Building validation windows from %s", valid_regions)
        x_valid, y_valid, _ = build_dataset(
            fasta_path, bigwig_path, valid_regions, seq_len=seq_len, step=step
        )
    else:
        logger.info(
            "No --valid_regions given; holding out %.0f%% of training windows",
            100 * valid_fraction,
        )
        x_train, y_train, x_valid, y_valid = _split_train_valid(
            x_train, y_train, valid_fraction, seed
        )
    logger.info(
        "Training on %d windows, validating on %d windows",
        x_train.shape[0],
        x_valid.shape[0],
    )

    train_loader = get_dataloader(
        DNASeqDataset(x_train, y_train),
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        # Drop a size-1 trailing batch: BatchNorm1d needs >1 sample in train mode.
        drop_last=True,
    )
    valid_loader = get_dataloader(
        DNASeqDataset(x_valid, y_valid),
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
    )

    # --- model / optimiser --------------------------------------------------
    dev = _select_device(device)
    logger.info("Training on device: %s", dev)
    model = EditNet(
        seq_len=seq_len, n_filters=n_filters, kernel_size=kernel_size
    ).to(dev)
    # The model outputs a positive Poisson rate (lambda) via Softplus, so the
    # loss receives lambda directly (log_input=False). full=False drops the
    # data-only Stirling term, which does not affect the gradients.
    loss_fn = nn.PoissonNLLLoss(log_input=False, full=False)
    optimizer = torch.optim.Adam(
        model.parameters(), lr=lr, weight_decay=weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", patience=lr_patience, min_lr=min_lr
    )

    # --- training loop ------------------------------------------------------
    out_path = os.path.join(out_dir, f"{out_name}.pth")
    train_losses: list[float] = []
    valid_losses: list[float] = []
    best_valid = float("inf")
    best_epoch = -1

    for epoch in range(epochs):
        train_loss = _run_epoch(model, train_loader, loss_fn, dev, optimizer)
        valid_loss = _run_epoch(model, valid_loader, loss_fn, dev)
        scheduler.step(valid_loss)
        train_losses.append(train_loss)
        valid_losses.append(valid_loss)

        improved = valid_loss < best_valid
        if improved:
            best_valid = valid_loss
            best_epoch = epoch
            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "config": model.config,
                    "train_losses": train_losses,
                    "valid_losses": valid_losses,
                    "best_epoch": best_epoch,
                    "best_valid_loss": best_valid,
                },
                out_path,
            )
        logger.info(
            "epoch %3d/%d  train_loss=%.5f  valid_loss=%.5f  lr=%.2e%s",
            epoch + 1,
            epochs,
            train_loss,
            valid_loss,
            optimizer.param_groups[0]["lr"],
            "  *" if improved else "",
        )

    logger.info(
        "Done. Best valid_loss=%.5f at epoch %d. Checkpoint: %s",
        best_valid,
        best_epoch + 1,
        out_path,
    )
    return out_path
