import numpy as np
import pytest

from exocam_accelerate.phase_space import (
    PhaseGateConfig,
    extrapolate,
    fit_hyperbolic,
    fit_linear,
    fit_saturating,
)


def test_linear_recovers_intercept():
    N = np.linspace(-5, -1, 20)
    X = 200.0 + 3.0 * N
    fit = fit_linear(N, X)
    assert fit.predict(0.0) == pytest.approx(200.0)
    assert fit.corr == pytest.approx(1.0)


def test_saturating_recovers_curve():
    N = np.linspace(-20, -1, 30)
    X = 250.0 - 15.0 * (1.0 - np.exp(N / 6.0))
    fit = fit_saturating(N, X)
    assert fit.params[2] == pytest.approx(6.0, rel=1e-3)
    assert fit.predict(0.0) == pytest.approx(250.0, abs=1e-6)


def test_hyperbolic_recovers_stefan_law():
    h = np.linspace(20, 40, 25)
    N = 0.2 - 60.0 / h
    fit = fit_hyperbolic(N, h)
    a, b = fit.params
    assert a == pytest.approx(0.2)
    assert b == pytest.approx(-60.0)
    assert fit.predict(-1.0) == pytest.approx(50.0)
    # past the asymptote the thickness is infinite
    assert np.isinf(fit.predict(0.5))


def test_hyperbolic_rejects_nonpositive():
    with pytest.raises(ValueError):
        fit_hyperbolic([-2, -1.5, -1], [1.0, 0.0, 2.0])


class TestExtrapolateGate:
    def _stefan(self, n=20, noise=0.0, seed=0):
        rng = np.random.default_rng(seed)
        h = np.linspace(30, 40, n)
        N = 0.1 - 70.0 / h + noise * rng.standard_normal(n)
        return N, h

    def test_accepts_clean_stefan(self):
        N, h = self._stefan()
        r = extrapolate("hi", N, h, -1.2, 40.0, "hyperbolic")
        assert r.accepted, r.reasons
        assert r.prediction == pytest.approx(70.0 / 1.3)

    def test_refuses_past_asymptote(self):
        N, h = self._stefan()
        r = extrapolate("hi", N, h, 0.2, 40.0, "hyperbolic",
                        PhaseGateConfig(max_extrapolation_ratio=100))
        assert not r.accepted
        assert any("asymptote" in s for s in r.reasons)
        assert np.isnan(r.prediction)

    def test_refuses_far_extrapolation(self):
        N, h = self._stefan()
        r = extrapolate("hi", N, h, -0.1, 40.0, "hyperbolic",
                        PhaseGateConfig(max_extrapolation_ratio=2.0))
        assert not r.accepted
        assert any("beyond the window" in s for s in r.reasons)

    def test_refuses_uncorrelated(self):
        N, h = self._stefan(noise=2.0)
        r = extrapolate("hi", N, h, -1.6, 40.0, "hyperbolic")
        assert not r.accepted
        assert any("conduction-limited" in s for s in r.reasons)

    def test_too_few_points(self):
        r = extrapolate("TS", [-2, -1.5], [210, 209], 0.0, 209.0)
        assert not r.accepted

    def test_temperature_kink_refused(self):
        # strongly curved: linear and saturating endpoints disagree
        N = np.linspace(-30, -10, 20)
        X = 250.0 - 40.0 * (1.0 - np.exp(N / 3.0))
        r = extrapolate("TS", N, X, 0.0, float(X[-1]), "linear",
                        PhaseGateConfig(max_extrapolation_ratio=10,
                                        max_relative_disagreement=0.2))
        assert not r.accepted
        assert any("disagree" in s for s in r.reasons)

    def test_temperature_linear_accepted(self):
        N = np.linspace(-3, -2, 20)
        X = 210.0 + 1.5 * N
        r = extrapolate("TS", N, X, -1.0, float(X[-1]), "linear")
        assert r.accepted, r.reasons
        assert r.prediction == pytest.approx(208.5)

    def test_unknown_form(self):
        with pytest.raises(ValueError):
            extrapolate("x", [1, 2], [1, 2], 0, 1, "cubic")
