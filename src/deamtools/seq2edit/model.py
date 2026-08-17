"""The seq2edit CNN model (``EditNet``) and its sequence/signal encoding.

``EditNet`` learns a *DNA sequence -> per-base editing* map (the deaminase
sequence bias): given a fixed-length window of one-hot DNA it predicts the
editing signal at every base of that window. Trained on naked/deproteinised DNA
(or any background where editing reflects enzyme preference rather than
chromatin state), it yields an *expected* track that footprint and occupancy
analyses divide out.

The architecture follows the ACCESS-ATAC ``cnn_bias_model``
(https://github.com/pinellolab/ACCESS-ATAC/tree/main/cnn_bias_model): two
``Conv1d -> ReLU -> MaxPool -> Dropout`` blocks followed by a fully-connected
head. Two deliberate changes:

* the flattened feature size feeding the first ``Linear`` is computed from the
  architecture (:func:`flatten_dim`) rather than hard-coded to ``928``, so the
  network works for window/filter sizes other than ``128 / 32 / 5``; and
* the head ends in a ``Softplus`` so the output is a strictly positive rate
  ``lambda >= 0`` -- the mean of a Poisson distribution over per-base edit
  counts. Training uses a Poisson negative-log-likelihood loss (see
  :mod:`deamtools.seq2edit.train`) rather than MSE, which matches the
  count-like nature of editing signal.

The encoding/data-prep helpers (:func:`one_hot_encode`, :func:`build_dataset`,
:func:`flatten_dim`, ...) do not use PyTorch at call time, but importing this
module pulls in ``torch`` for the model. Install the extra with
``pip install 'deamtools[seq2edit]'``.
"""

from __future__ import annotations

import logging

import numpy as np
import pyBigWig
import pysam
from torch import nn

logger = logging.getLogger(__name__)

# Column order of the one-hot encoding.
BASES: tuple[str, ...] = ("A", "C", "G", "T")
_BASE_TO_COL = {b: i for i, b in enumerate(BASES)}


# --------------------------------------------------------------------------- #
# Sequence / signal encoding.
# --------------------------------------------------------------------------- #
def one_hot_encode(seq: str) -> np.ndarray:
    """One-hot encode a DNA sequence as an ``(len(seq), 4)`` float32 array.

    Columns are ordered ``A, C, G, T`` (see :data:`BASES`). The input is
    upper-cased first; any base that is not ``A/C/G/T`` (e.g. ``N`` or a soft-
    masked lowercase letter that is not a canonical base) is encoded as an
    all-zero row, contributing no signal rather than a spurious base. This is
    the one deliberate deviation from the ACCESS-ATAC reference, which maps
    ``N`` to ``T``.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(len(seq), 4)`` and dtype ``float32``.
    """
    vec = np.zeros((len(seq), 4), dtype=np.float32)
    for i, base in enumerate(seq.upper()):
        col = _BASE_TO_COL.get(base)
        if col is not None:
            vec[i, col] = 1.0
    return vec


def conv1d_out_len(length: int, kernel_size: int, pool: int = 2) -> int:
    """Length after one ``Conv1d`` (stride 1, no padding) then ``MaxPool1d``.

    ``Conv1d`` with ``kernel_size`` and stride 1 maps ``L`` to
    ``L - (kernel_size - 1)``; ``MaxPool1d(pool)`` then floor-divides by
    ``pool``. Mirrors the two convolutional blocks of :class:`EditNet`.
    """
    return (length - (kernel_size - 1)) // pool


def flatten_dim(seq_len: int, n_filters: int, kernel_size: int) -> int:
    """Flattened feature size after the two conv blocks of :class:`EditNet`.

    Computed from the architecture rather than hard-coded (the ACCESS-ATAC
    reference hard-codes ``928`` for ``seq_len=128, n_filters=32,
    kernel_size=5``) so the network generalises to other window/filter sizes.
    """
    after1 = conv1d_out_len(seq_len, kernel_size)
    after2 = conv1d_out_len(after1, kernel_size)
    if after2 <= 0:
        raise ValueError(
            f"seq_len={seq_len} is too short for kernel_size={kernel_size}: "
            "the second convolution leaves no positions. Increase --seq_len "
            "or decrease --kernel_size."
        )
    return n_filters * after2


def iter_windows(
    chrom: str, start: int, end: int, seq_len: int, step: int
) -> list[tuple[int, int]]:
    """Tile ``[start, end)`` into ``seq_len``-wide windows spaced by ``step``.

    A trailing stretch shorter than ``seq_len`` is dropped so every window has
    exactly ``seq_len`` bases. Returns a list of ``(win_start, win_end)``.
    """
    windows: list[tuple[int, int]] = []
    pos = start
    while pos + seq_len <= end:
        windows.append((pos, pos + seq_len))
        pos += step
    return windows


def build_dataset(
    fasta_path: str,
    bigwig_path: str,
    bed_path: str,
    seq_len: int = 128,
    step: int | None = None,
) -> tuple[np.ndarray, np.ndarray, list[tuple[str, int, int]]]:
    """Build training arrays by tiling BED regions into fixed windows.

    For every ``seq_len``-wide window tiled across the BED intervals, the
    reference sequence is one-hot encoded (input ``x``) and the per-base
    BigWig signal over the same coordinates is read as the target (``y``).
    ``NaN`` BigWig values (uncovered bases) are treated as ``0``.

    Parameters
    ----------
    fasta_path : str
        Reference FASTA, indexed with ``samtools faidx`` (``.fai`` required).
    bigwig_path : str
        Per-base editing-signal BigWig (e.g. from ``deamtools bam2bw``),
        ideally from a naked/deproteinised-DNA control so the model captures
        enzyme sequence bias rather than chromatin state. The targets are
        treated as Poisson-distributed counts during training.
    bed_path : str
        BED file of regions to draw windows from.
    seq_len : int, optional
        Window width in bp (model input/output length). Default ``128``.
    step : int, optional
        Spacing between consecutive window starts. Defaults to ``seq_len``
        (non-overlapping tiling). Use a smaller value to oversample with
        overlapping windows.

    Returns
    -------
    x : numpy.ndarray
        One-hot inputs, shape ``(n_windows, seq_len, 4)``, dtype float32.
    y : numpy.ndarray
        Per-base targets, shape ``(n_windows, seq_len)``, dtype float32.
    coords : list of (str, int, int)
        ``(chrom, start, end)`` for each window, in row order.
    """
    if step is None:
        step = seq_len
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")

    # Local import to reuse the project's merged-interval BED loader.
    from deamtools.utils import _load_regions

    regions = _load_regions(bed_path)

    fasta = pysam.FastaFile(fasta_path)
    bw = pyBigWig.open(bigwig_path)
    try:
        bw_chroms = set(bw.chroms().keys())
        fa_chroms = set(fasta.references)

        x_list: list[np.ndarray] = []
        y_list: list[np.ndarray] = []
        coords: list[tuple[str, int, int]] = []
        skipped = 0

        for row in regions.itertuples(index=False):
            chrom, r_start, r_end = row.chrom, int(row.start), int(row.end)
            if chrom not in fa_chroms or chrom not in bw_chroms:
                skipped += 1
                continue
            chrom_len = fasta.get_reference_length(chrom)
            for w_start, w_end in iter_windows(
                chrom, r_start, min(r_end, chrom_len), seq_len, step
            ):
                seq = fasta.fetch(chrom, w_start, w_end)
                if len(seq) != seq_len:
                    skipped += 1
                    continue
                signal = bw.values(chrom, w_start, w_end, numpy=True)
                signal = np.nan_to_num(signal, nan=0.0).astype(np.float32)
                x_list.append(one_hot_encode(seq))
                y_list.append(signal)
                coords.append((chrom, w_start, w_end))
    finally:
        bw.close()
        fasta.close()

    if not x_list:
        raise ValueError(
            f"No usable {seq_len}-bp windows found in {bed_path}. Check that "
            "the BED regions are at least seq_len wide and that chromosome "
            "names match between the FASTA, BigWig, and BED."
        )

    x = np.stack(x_list).astype(np.float32)
    y = np.stack(y_list).astype(np.float32)
    logger.info(
        "Built %d windows of %d bp (%d region(s) skipped) from %s",
        x.shape[0],
        seq_len,
        skipped,
        bed_path,
    )
    return x, y, coords


# --------------------------------------------------------------------------- #
# The network.
# --------------------------------------------------------------------------- #
class EditNet(nn.Module):
    """CNN mapping a one-hot DNA window to a per-base editing rate.

    The output is a strictly positive rate ``lambda`` (via ``Softplus``) — the
    mean of a Poisson over per-base edit counts — so it pairs with the Poisson
    negative-log-likelihood loss used for training.

    Parameters
    ----------
    seq_len : int
        Input/output window length in bp. The network consumes a
        ``(batch, seq_len, 4)`` one-hot tensor and predicts ``(batch, seq_len)``.
    n_filters : int
        Number of convolutional filters in both conv blocks.
    kernel_size : int
        Convolution kernel width.
    hidden : int
        Width of the two fully-connected hidden layers.
    dropout : float
        Dropout probability used throughout.
    """

    def __init__(
        self,
        seq_len: int = 128,
        n_filters: int = 32,
        kernel_size: int = 5,
        hidden: int = 1024,
        dropout: float = 0.25,
    ) -> None:
        super().__init__()
        self.seq_len = seq_len
        self.n_filters = n_filters
        self.kernel_size = kernel_size

        self.conv1 = nn.Sequential(
            nn.Conv1d(4, n_filters, kernel_size=kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(n_filters, n_filters, kernel_size=kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),
        )

        flat = flatten_dim(seq_len, n_filters, kernel_size)
        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(flat, hidden),
            nn.ReLU(),
            nn.BatchNorm1d(hidden),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.BatchNorm1d(hidden),
            nn.Dropout(dropout),
            nn.Linear(hidden, seq_len),
        )
        # Map to a strictly positive Poisson rate lambda for every base.
        self.output = nn.Softplus()

    def forward(self, x):  # noqa: D102 - standard nn.Module forward
        # x: (batch, seq_len, 4) -> (batch, 4, seq_len) for Conv1d.
        x = x.permute(0, 2, 1)
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.fc(x)
        return self.output(x)

    @property
    def config(self) -> dict:
        """Hyper-parameters needed to rebuild this model for inference."""
        return {
            "seq_len": self.seq_len,
            "n_filters": self.n_filters,
            "kernel_size": self.kernel_size,
        }
