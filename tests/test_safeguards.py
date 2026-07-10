import numpy as np
import pytest

from exocam_accelerate.safeguards import (
    GateConfig,
    assess_trustworthiness,
    clip_step,
)
from exocam_accelerate.trends import TrendSeries


def linear_series(slope=-0.5, n=24, noise=0.0, seed=0):
    t = np.linspace(0, 10, n)
    rng = np.random.default_rng(seed)
    return TrendSeries(times=t, values=250.0 + slope * t + rng.normal(0, noise, n))


class TestGateAccepts:
    def test_clean_linear_drift_passes(self):
        result = assess_trustworthiness(linear_series())
        assert result.ok
        assert bool(result) is True
        assert result.reasons == ()

    def test_noisy_linear_drift_passes_with_tolerance(self):
        # small noise atop a strong drift: curvature contribution stays small
        result = assess_trustworthiness(
            linear_series(slope=-1.0, noise=0.05, seed=3),
            GateConfig(negligible_change=0.5),
        )
        assert result.ok

    def test_converged_flat_series_passes_as_converged(self):
        t = np.linspace(0, 10, 24)
        rng = np.random.default_rng(1)
        s = TrendSeries(times=t, values=280.0 + rng.normal(0, 1e-4, len(t)))
        result = assess_trustworthiness(s, GateConfig(negligible_change=0.01))
        assert result.ok
        assert result.converged_mask.all()


class TestGateRefuses:
    def test_strong_curvature_refuses(self):
        # accelerating drift — approach to a feedback threshold
        t = np.linspace(0, 10, 24)
        s = TrendSeries(times=t, values=250.0 - 0.1 * t - 0.2 * t**2)
        result = assess_trustworthiness(s)
        assert not result.ok
        assert any("curvature" in r for r in result.reasons)

    def test_sign_flip_refuses(self):
        # V-shaped: system reversed direction inside the window
        t = np.linspace(0, 10, 24)
        s = TrendSeries(times=t, values=250.0 + np.abs(t - 5.0))
        result = assess_trustworthiness(s)
        assert not result.ok
        assert any("sign" in r for r in result.reasons)

    def test_pure_curvature_on_flat_trend_refuses_by_default(self):
        # symmetric parabola: zero net linear change but strong curvature;
        # with the conservative default negligible_change=0 this must refuse
        t = np.linspace(0, 10, 24)
        s = TrendSeries(times=t, values=250.0 + 0.5 * (t - 5.0) ** 2)
        result = assess_trustworthiness(s)
        assert not result.ok

    def test_per_layer_mixed_verdict_refuses_overall(self):
        t = np.linspace(0, 10, 24)
        good = 250.0 - 0.5 * t
        bad = 250.0 - 0.1 * t - 0.2 * t**2
        s = TrendSeries(times=t, values=np.stack([good, bad], axis=1))
        result = assess_trustworthiness(s)
        assert not result.ok
        assert result.ok_mask[0] and not result.ok_mask[1]


class TestClip:
    def test_below_limit_untouched(self):
        delta = np.array([1.0, -2.0, 0.5])
        report = clip_step(delta, 5.0)
        assert np.array_equal(report.clipped, delta)
        assert not report.any_clipped

    def test_clips_and_preserves_sign(self):
        delta = np.array([80.0, -120.0, 10.0])
        report = clip_step(delta, 50.0)   # Turbet's 50 K
        assert np.array_equal(report.clipped, [50.0, -50.0, 10.0])
        assert report.n_clipped == 2
        assert report.max_requested == 120.0
        assert report.limit == 50.0

    def test_nonpositive_limit_raises(self):
        with pytest.raises(ValueError):
            clip_step(np.array([1.0]), 0.0)
