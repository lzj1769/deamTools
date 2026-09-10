"""Tests for the ChromBPNet-architecture model in ``deamtools.seq2edit.bpnet``."""

from __future__ import annotations

import math

import pytest

torch = pytest.importorskip("torch")

from deamtools.seq2edit.bpnet import (  # noqa: E402
    BPNet,
    BPNetLoss,
    body_out_len,
    count_loss,
    estimate_counts_loss_weight,
    multinomial_nll,
    profile_out_len,
    required_input_len,
)


class TestLengthArithmetic:
    def test_chrombpnet_defaults_give_2114_to_1000(self):
        """The published ChromBPNet configuration must reproduce exactly."""
        assert profile_out_len(2114) == 1000
        assert required_input_len(1000) == 2114

    def test_body_matches_step_by_step_trim(self):
        length = 2114 - (21 - 1)
        for i in range(1, 9):
            length -= 2 * (2**i)
        assert body_out_len(2114, conv1_kernel_size=21, n_dil_layers=8) == length

    def test_required_input_len_inverts_profile_out_len(self):
        for out_len in (100, 250, 1000, 1024):
            for n_dil in (4, 6, 8):
                n = required_input_len(out_len, n_dil_layers=n_dil)
                assert profile_out_len(n, n_dil_layers=n_dil) == out_len

    def test_dilated_stack_trim_closed_form(self):
        """`2**(n+2) - 4` must equal the summed per-layer trim."""
        for n in range(1, 11):
            assert sum(2 * 2**i for i in range(1, n + 1)) == 2 ** (n + 2) - 4


class TestConstruction:
    def test_rejects_too_short_input(self):
        with pytest.raises(ValueError, match="short of output_len"):
            BPNet(input_len=500, output_len=1000, filters=4)

    def test_error_names_a_usable_input_len(self):
        with pytest.raises(ValueError) as exc:
            BPNet(input_len=500, output_len=1000, filters=4)
        assert "input_len=2114" in str(exc.value)

    def test_rejects_odd_crop(self):
        """An odd surplus cannot be centred, so it is refused up front."""
        base = required_input_len(100, n_dil_layers=2, profile_kernel_size=5)
        with pytest.raises(ValueError, match="odd difference"):
            BPNet(
                input_len=base + 1,
                output_len=100,
                filters=4,
                n_dil_layers=2,
                conv1_kernel_size=3,
                profile_kernel_size=5,
            )

    def test_even_surplus_is_accepted(self):
        base = required_input_len(100, n_dil_layers=2, profile_kernel_size=5)
        model = BPNet(
            input_len=base + 2,
            output_len=100,
            filters=4,
            n_dil_layers=2,
            conv1_kernel_size=3,
            profile_kernel_size=5,
        )
        assert model.output_len == 100

    def test_dilations_double(self):
        model = BPNet(
            input_len=required_input_len(20, n_dil_layers=4, profile_kernel_size=5),
            output_len=20,
            filters=4,
            n_dil_layers=4,
            profile_kernel_size=5,
        )
        assert [c.dilation[0] for c in model.dilated] == [2, 4, 8, 16]


def _small_model(**kw):
    """A tiny model with the same topology, for fast tests."""
    params = dict(
        output_len=20,
        filters=8,
        n_dil_layers=3,
        conv1_kernel_size=5,
        profile_kernel_size=7,
        n_tasks=1,
    )
    params.update(kw)
    params["input_len"] = required_input_len(
        params["output_len"],
        params["conv1_kernel_size"],
        params["n_dil_layers"],
        params["profile_kernel_size"],
    )
    return BPNet(**params)


class TestForward:
    def test_output_shapes(self):
        model = _small_model()
        x = torch.rand(3, model.input_len, 4)
        profile, log_counts = model(x)
        assert profile.shape == (3, 1, 20)
        assert log_counts.shape == (3, 1)

    def test_multi_task_shapes(self):
        model = _small_model(n_tasks=4)
        x = torch.rand(2, model.input_len, 4)
        profile, log_counts = model(x)
        assert profile.shape == (2, 4, 20)
        assert log_counts.shape == (2, 4)

    def test_profile_is_logits_not_a_distribution(self):
        """The head must emit raw logits; softmax is the caller's job."""
        model = _small_model()
        x = torch.rand(4, model.input_len, 4)
        profile, _ = model(x)
        sums = profile.softmax(dim=-1).sum(dim=-1)
        assert torch.allclose(sums, torch.ones_like(sums), atol=1e-5)
        assert not torch.allclose(profile.sum(dim=-1), torch.ones_like(sums), atol=1e-2)

    def test_predict_signal_totals_match_counts_head(self):
        model = _small_model()
        x = torch.rand(3, model.input_len, 4)
        _, log_counts = model(x)
        signal = model.predict_signal(x)
        expected = torch.expm1(log_counts).clamp_min(0.0)
        assert torch.allclose(signal.sum(dim=-1), expected, atol=1e-4)

    def test_gradients_reach_every_layer(self):
        model = _small_model()
        x = torch.rand(2, model.input_len, 4)
        profile, log_counts = model(x)
        (profile.sum() + log_counts.sum()).backward()
        missing = [n for n, p in model.named_parameters() if p.grad is None]
        assert missing == []

    def test_config_round_trips(self):
        model = _small_model(n_tasks=2)
        rebuilt = BPNet(**model.config)
        assert rebuilt.config == model.config
        rebuilt.load_state_dict(model.state_dict())


class TestLosses:
    def test_multinomial_nll_is_minimised_by_the_true_shape(self):
        counts = torch.tensor([[[0.0, 10.0, 30.0, 60.0]]])
        perfect = torch.log(torch.tensor([[[1e-8, 0.1, 0.3, 0.6]]]))
        wrong = torch.log(torch.tensor([[[0.6, 0.3, 0.1, 1e-8]]]))
        assert multinomial_nll(perfect, counts) < multinomial_nll(wrong, counts)

    def test_multinomial_nll_ignores_overall_scale(self):
        """Shifting logits by a constant is a no-op after softmax."""
        counts = torch.tensor([[[1.0, 2.0, 7.0]]])
        logits = torch.tensor([[[0.5, -0.2, 1.3]]])
        assert torch.allclose(
            multinomial_nll(logits, counts), multinomial_nll(logits + 4.0, counts)
        )

    def test_multinomial_nll_matches_hand_computation(self):
        counts = torch.tensor([[2.0, 3.0]])
        logits = torch.tensor([[0.0, 0.0]])  # p = 0.5, 0.5
        # log C(5;2,3) + 5*log(0.5)
        expected = -(math.log(math.comb(5, 2)) + 5 * math.log(0.5))
        assert multinomial_nll(logits, counts).item() == pytest.approx(
            expected, rel=1e-5
        )

    def test_count_loss_is_zero_when_exact(self):
        counts = torch.tensor([[[1.0, 2.0, 3.0]]])
        log_counts = torch.log1p(counts.sum(dim=-1))
        assert count_loss(log_counts, counts).item() == pytest.approx(0.0, abs=1e-6)

    def test_combined_loss_reports_components(self):
        model = _small_model()
        x = torch.rand(4, model.input_len, 4)
        counts = torch.randint(0, 5, (4, 1, 20)).float()
        profile, log_counts = model(x)
        total, prof_term, count_term = BPNetLoss(counts_weight=3.0)(
            profile, log_counts, counts
        )
        assert torch.allclose(total, prof_term + 3.0 * count_term)

    def test_counts_weight_scales_only_the_counts_term(self):
        model = _small_model()
        x = torch.rand(4, model.input_len, 4)
        counts = torch.randint(0, 5, (4, 1, 20)).float()
        profile, log_counts = model(x)
        _, p1, c1 = BPNetLoss(counts_weight=1.0)(profile, log_counts, counts)
        _, p2, c2 = BPNetLoss(counts_weight=9.0)(profile, log_counts, counts)
        assert torch.allclose(p1, p2)
        assert torch.allclose(c1, c2)

    def test_estimate_counts_loss_weight_uses_median_total(self):
        counts = torch.zeros(5, 1, 4)
        counts[:, 0, 0] = torch.tensor([10.0, 20.0, 30.0, 40.0, 50.0])
        assert estimate_counts_loss_weight(counts, scale=10.0) == pytest.approx(3.0)


class TestTrainability:
    def test_overfits_a_single_window(self):
        """A few steps on one example must reduce the loss substantially."""
        torch.manual_seed(0)
        model = _small_model()
        x = torch.rand(2, model.input_len, 4)
        counts = torch.zeros(2, 1, 20)
        counts[:, 0, 5:10] = 20.0  # a sharp, learnable peak

        loss_fn = BPNetLoss(counts_weight=1.0)
        opt = torch.optim.Adam(model.parameters(), lr=1e-2)

        profile, log_counts = model(x)
        first = loss_fn(profile, log_counts, counts)[0].item()
        for _ in range(60):
            opt.zero_grad()
            profile, log_counts = model(x)
            loss = loss_fn(profile, log_counts, counts)[0]
            loss.backward()
            opt.step()
        assert loss.item() < first * 0.5
