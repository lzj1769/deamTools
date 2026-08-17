"""Tests for deamtools.seq2edit.

``deamtools.seq2edit.model`` imports PyTorch (for ``EditNet``), so the whole
module — including the otherwise torch-free encoding/sizing helpers it also
hosts — is skipped when the optional ``torch`` extra is not installed.
"""

import os

import numpy as np
import pyBigWig
import pysam
import pytest

pytest.importorskip("torch")

from deamtools.seq2edit.model import (  # noqa: E402 - after importorskip guard
    build_dataset,
    conv1d_out_len,
    flatten_dim,
    iter_windows,
    one_hot_encode,
)


# --------------------------------------------------------------------------- #
# Fixtures: a tiny synthetic FASTA + BigWig.
# --------------------------------------------------------------------------- #
def _write_fasta(path: str, seqs: dict[str, str]) -> None:
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")
    pysam.faidx(path)


def _write_bigwig(path: str, chrom_len: int, values: list[float]) -> None:
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", chrom_len)])
    positions = list(range(len(values)))
    bw.addEntries(
        ["chr1"] * len(values),
        positions,
        ends=[p + 1 for p in positions],
        values=[float(v) for v in values],
    )
    bw.close()


def _write_bed(path: str, rows: list[tuple[str, int, int]]) -> None:
    with open(path, "w") as f:
        for chrom, start, end in rows:
            f.write(f"{chrom}\t{start}\t{end}\n")


class TestOneHotEncode:
    def test_canonical_bases(self):
        x = one_hot_encode("ACGT")
        expected = np.eye(4, dtype=np.float32)
        np.testing.assert_array_equal(x, expected)

    def test_lowercase_is_upcased(self):
        np.testing.assert_array_equal(one_hot_encode("acgt"), one_hot_encode("ACGT"))

    def test_n_is_all_zero(self):
        x = one_hot_encode("ANG")
        np.testing.assert_array_equal(x[1], np.zeros(4, dtype=np.float32))
        assert x.dtype == np.float32
        assert x.shape == (3, 4)


class TestSizing:
    def test_conv1d_out_len(self):
        # L=128, kernel 5: (128-4)//2 = 62
        assert conv1d_out_len(128, 5) == 62
        # then (62-4)//2 = 29
        assert conv1d_out_len(62, 5) == 29

    def test_flatten_dim_matches_reference(self):
        # ACCESS-ATAC hard-codes 928 for seq_len=128, n_filters=32, kernel=5.
        assert flatten_dim(128, 32, 5) == 928

    def test_flatten_dim_too_short_raises(self):
        with pytest.raises(ValueError):
            flatten_dim(8, 32, 5)


class TestIterWindows:
    def test_non_overlapping_tiling(self):
        wins = iter_windows("chr1", 0, 256, seq_len=128, step=128)
        assert wins == [(0, 128), (128, 256)]

    def test_trailing_short_window_dropped(self):
        wins = iter_windows("chr1", 0, 200, seq_len=128, step=128)
        assert wins == [(0, 128)]

    def test_overlapping_step(self):
        wins = iter_windows("chr1", 0, 200, seq_len=128, step=64)
        assert wins == [(0, 128), (64, 192)]


class TestBuildDataset:
    def test_shapes_and_targets(self, tmp_path):
        seq_len = 16
        chrom_len = 64
        rng = np.random.default_rng(0)
        seq = "".join(rng.choice(list("ACGT"), size=chrom_len))
        signal = list(np.arange(chrom_len, dtype=float))  # 0,1,2,... per base

        fasta = str(tmp_path / "ref.fa")
        bw = str(tmp_path / "sig.bw")
        bed = str(tmp_path / "regions.bed")
        _write_fasta(fasta, {"chr1": seq})
        _write_bigwig(bw, chrom_len, signal)
        _write_bed(bed, [("chr1", 0, 32)])  # -> two 16-bp windows

        x, y, coords = build_dataset(fasta, bw, bed, seq_len=seq_len)

        assert x.shape == (2, seq_len, 4)
        assert y.shape == (2, seq_len)
        assert coords == [("chr1", 0, 16), ("chr1", 16, 32)]
        # Target equals the per-base BigWig signal over the window.
        np.testing.assert_array_almost_equal(y[0], np.arange(0, 16))
        np.testing.assert_array_almost_equal(y[1], np.arange(16, 32))
        # Input one-hot round-trips to the reference sequence.
        bases = np.array(list("ACGT"))
        decoded = "".join(bases[x[0].argmax(axis=1)])
        assert decoded == seq[:16]

    def test_no_windows_raises(self, tmp_path):
        fasta = str(tmp_path / "ref.fa")
        bw = str(tmp_path / "sig.bw")
        bed = str(tmp_path / "regions.bed")
        _write_fasta(fasta, {"chr1": "ACGT" * 8})
        _write_bigwig(bw, 32, [0.0] * 32)
        _write_bed(bed, [("chr1", 0, 10)])  # shorter than seq_len=128
        with pytest.raises(ValueError):
            build_dataset(fasta, bw, bed, seq_len=128)


class TestModelAndTrain:
    def test_editnet_forward_shape(self):
        import torch

        from deamtools.seq2edit.model import EditNet

        model = EditNet(seq_len=128)
        model.eval()  # disable dropout/batchnorm randomness for the assertion
        x = torch.zeros((4, 128, 4))
        out = model(x)
        assert out.shape == (4, 128)
        # Softplus head emits a strictly positive Poisson rate.
        assert (out > 0).all()

    def test_run_train_writes_checkpoint(self, tmp_path):
        import torch

        from deamtools.seq2edit.train import run_train

        seq_len = 16
        chrom_len = 512
        rng = np.random.default_rng(1)
        seq = "".join(rng.choice(list("ACGT"), size=chrom_len))
        signal = list(rng.random(chrom_len))

        fasta = str(tmp_path / "ref.fa")
        bw = str(tmp_path / "sig.bw")
        bed = str(tmp_path / "regions.bed")
        _write_fasta(fasta, {"chr1": seq})
        _write_bigwig(bw, chrom_len, signal)
        _write_bed(bed, [("chr1", 0, chrom_len)])  # 32 windows of 16 bp

        out_path = run_train(
            fasta_path=fasta,
            bigwig_path=bw,
            train_regions=bed,
            out_dir=str(tmp_path / "model"),
            out_name="bias",
            seq_len=seq_len,
            epochs=2,
            batch_size=8,
            device="cpu",
            seed=0,
        )

        assert os.path.exists(out_path)
        ckpt = torch.load(out_path, map_location="cpu", weights_only=False)
        assert ckpt["config"]["seq_len"] == seq_len
        assert len(ckpt["train_losses"]) == 2
        assert "model_state_dict" in ckpt
