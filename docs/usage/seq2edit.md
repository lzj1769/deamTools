# seq2edit

Model the **deaminase sequence bias** by learning a map from DNA sequence to
per-base editing signal. A deaminase does not edit every accessible cytosine
equally — it has a sequence preference (e.g. DddA's `TC` context). `seq2edit`
trains a convolutional neural network to predict, for a window of one-hot DNA,
the editing each base *would* receive on the basis of sequence alone. Trained on
a naked/deproteinised-DNA control, this gives an **expected** track that
downstream footprint and occupancy analyses can divide out to separate enzyme
bias from genuine protein protection.

The design follows the ACCESS-ATAC
[`cnn_bias_model`](https://github.com/pinellolab/ACCESS-ATAC/tree/main/cnn_bias_model).

The module has three stages — **train** (implemented here), **predict**, and
**interpret** (planned). This page documents `seq2edit train`.

:::{note}
`seq2edit` needs the optional **PyTorch** dependency, which is not installed by
default:

```bash
pip install 'deamtools[seq2edit]'
```

The data-preparation helpers work without it; only model training/inference
require `torch`.
:::

## Synopsis

```
deamtools seq2edit train --fasta FILE --bigwig FILE --train_regions FILE \
    --out_dir DIR --out_name NAME [options]
```

## Required inputs

| Argument | Description |
|---|---|
| `--fasta FILE` | Reference FASTA indexed with `samtools faidx` (`.fai` required). Source of the input sequence. |
| `--bigwig FILE` | Per-base editing-signal BigWig (the regression target), e.g. from [`deamtools bam2bw`](bam2bw.md). Use a **naked/deproteinised-DNA** control so the model learns enzyme bias rather than chromatin state. |
| `--train_regions FILE` | BED of regions tiled into training windows. |
| `--out_dir DIR` | Output directory. Created automatically if it does not exist. |
| `--out_name NAME` | Base name for the checkpoint. Writes `<out_dir>/<out_name>.pth`. |

## Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--valid_regions FILE` | *(hold-out split)* | BED of validation regions. If omitted, a random `--valid_fraction` of the training windows is held out. |
| `--valid_fraction FLOAT` | `0.1` | Hold-out fraction used when `--valid_regions` is not given. |
| `--seq_len INT` | `128` | Window width in bp — the model's input and output length. |
| `--step INT` | `--seq_len` | Spacing between window starts. Use a smaller value to oversample with overlapping windows. |
| `--n_filters INT` | `32` | Convolutional filters per conv block. |
| `--kernel_size INT` | `5` | Convolution kernel width. |
| `--epochs INT` | `200` | Number of training epochs. |
| `--batch_size INT` | `512` | Mini-batch size. |
| `--lr FLOAT` | `3e-4` | Adam learning rate. |
| `--weight_decay FLOAT` | `1e-4` | Adam weight decay (L2). |
| `--num_workers INT` | `0` | DataLoader worker processes. |
| `--seed INT` | `42` | RNG seed for reproducibility. |
| `--device STR` | *(auto)* | Compute device (`cpu`, `cuda`, `mps`). Auto-detects CUDA, then Apple MPS, then CPU. |
| `--log_level LEVEL` | `INFO` | Global flag (before the subcommand): `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

## How it works

**Windowing.** Each interval in `--train_regions` is tiled into non-overlapping
`--seq_len`-wide windows (a trailing stretch shorter than `seq_len` is dropped).
For every window:

- the reference sequence is one-hot encoded as the input `x` of shape
  `(seq_len, 4)` with column order `A, C, G, T` (non-`ACGT` bases such as `N`
  become an all-zero row); and
- the per-base `--bigwig` signal over the same coordinates is the target `y` of
  shape `(seq_len,)` (uncovered `NaN` bases are read as `0`).

**Model (`EditNet`).** Two `Conv1d → ReLU → MaxPool → Dropout` blocks feed a
fully-connected head that regresses the per-base editing vector for the window,
with a final `Softplus` so the output is a strictly positive rate `λ` — the mean
of a Poisson over per-base edit counts:

```
one-hot (seq_len, 4)
   → conv block 1 (n_filters, kernel_size)
   → conv block 2 (n_filters, kernel_size)
   → flatten → FC(1024) → FC(1024) → FC(seq_len) → Softplus
   → predicted per-base editing rate λ (seq_len,)
```

The flattened size feeding the first dense layer is derived from the
architecture, so non-default `--seq_len`/`--n_filters`/`--kernel_size` work
without code changes (the upstream reference hard-codes it for the `128 / 32 / 5`
defaults).

**Optimisation.** Editing signal is count-like, so the model is fit with a
**Poisson** negative-log-likelihood loss (`PoissonNLLLoss`, on the predicted
rate `λ`) rather than MSE — it maximises the Poisson likelihood of the observed
per-base counts. Adam (`--lr` `3e-4`, `--weight_decay` `1e-4`) with a
`ReduceLROnPlateau` schedule (patience 10, floor `1e-5`) on the validation loss.
Each epoch reports train/validation loss and the current learning rate; the
checkpoint is rewritten whenever the validation loss improves.

## Outputs

| File | Description |
|---|---|
| `<out_dir>/<out_name>.pth` | Best-validation checkpoint: model weights (`model_state_dict`), the `config` needed to rebuild the network (`seq_len`, `n_filters`, `kernel_size`), the per-epoch `train_losses`/`valid_losses`, and `best_epoch`/`best_valid_loss`. |

## Examples

```bash
# Train on naked-DNA editing signal over a set of regions
deamtools seq2edit train \
    --fasta hg38.fa --bigwig naked.bw \
    --train_regions regions.bed \
    --out_dir models --out_name bias

# Explicit validation regions, larger windows, fewer epochs, on GPU
deamtools seq2edit train \
    --fasta hg38.fa --bigwig naked.bw \
    --train_regions train.bed --valid_regions valid.bed \
    --seq_len 128 --epochs 50 --device cuda \
    --out_dir models --out_name bias
```

## Notes

- The `--bigwig` target should come from a control where editing reflects
  **sequence preference**, not chromatin state — typically naked or
  deproteinised DNA — otherwise the model conflates enzyme bias with
  accessibility.
- Windows are tiled, so making `--train_regions` broad (e.g. many peaks or a
  whole chromosome) gives the model more examples; `--step` smaller than
  `--seq_len` adds overlapping windows for additional coverage.
- The saved `config` lets later `predict`/`interpret` stages rebuild the exact
  network from the checkpoint.
