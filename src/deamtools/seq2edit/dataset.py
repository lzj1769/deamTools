"""PyTorch ``Dataset``/``DataLoader`` wrappers for seq2edit training.

Thin adapters over the numpy arrays produced by
:func:`deamtools.seq2edit.model.build_dataset`.
"""

from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset


class DNASeqDataset(Dataset):
    """Holds one-hot inputs ``x`` and optional per-base targets ``y``.

    Parameters
    ----------
    x : numpy.ndarray
        One-hot inputs, shape ``(n, seq_len, 4)``.
    y : numpy.ndarray, optional
        Per-base targets, shape ``(n, seq_len)``. When ``None`` (inference)
        ``__getitem__`` returns the input only.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray | None = None) -> None:
        self.x = torch.from_numpy(np.ascontiguousarray(x, dtype=np.float32))
        self.y = (
            None
            if y is None
            else torch.from_numpy(np.ascontiguousarray(y, dtype=np.float32))
        )

    def __len__(self) -> int:
        return self.x.shape[0]

    def __getitem__(self, idx: int):
        if self.y is None:
            return self.x[idx]
        return self.x[idx], self.y[idx]


def get_dataloader(
    dataset: DNASeqDataset,
    batch_size: int = 512,
    shuffle: bool = False,
    num_workers: int = 0,
    drop_last: bool = False,
) -> DataLoader:
    """Build a ``DataLoader`` with sensible defaults for this workload."""
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        drop_last=drop_last,
        pin_memory=torch.cuda.is_available(),
        persistent_workers=num_workers > 0,
    )
