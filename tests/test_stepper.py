import numpy as np
import pytest

from exocam_accelerate.safeguards import GateConfig
from exocam_accelerate.stepper import propose_step
from exocam_accelerate.trends import TrendSeries


def cooling_series(n=24, slope=-0.4):
    t = np.linspace(0, 10, n)
    return TrendSeries(times=t, values=260.0 + slope * t)


class TestProposeStep:
    def test_accepted_linear_step(self):
        s = cooling_series(slope=-0.4)
        p = propose_step(s, dt_years=50.0, max_abs_step=50.0)
        assert p.accepted and bool(p)
        assert np.isclose(p.tendency, -0.4)
        assert np.isclose(p.delta, -20.0)   # -0.4 K/yr * 50 yr
        assert not p.clip.any_clipped

    def test_clip_engages(self):
        s = cooling_series(slope=-2.0)
        p = propose_step(s, dt_years=100.0, max_abs_step=50.0)
        assert p.accepted
        assert np.isclose(p.delta, -50.0)   # requested -200, clipped
        assert p.clip.any_clipped

    def test_refusal_returns_no_delta(self):
        t = np.linspace(0, 10, 24)
        s = TrendSeries(times=t, values=260.0 - 0.1 * t - 0.3 * t**2)
        p = propose_step(s, dt_years=50.0, max_abs_step=50.0)
        assert not p.accepted and not bool(p)
        assert p.delta is None and p.tendency is None and p.clip is None
        assert p.gate.reasons

    def test_per_layer_delta(self):
        t = np.linspace(0, 10, 24)
        slopes = np.array([-0.5, -0.1, 0.2])
        values = 260.0 + t[:, None] * slopes[None, :]
        s = TrendSeries(times=t, values=values)
        p = propose_step(s, dt_years=10.0, max_abs_step=50.0)
        assert p.accepted
        assert np.allclose(p.delta, slopes * 10.0)

    def test_zero_dt_gives_zero_step(self):
        # Turbet's schedule ends at n=0: accelerate by nothing
        p = propose_step(cooling_series(), dt_years=0.0, max_abs_step=50.0)
        assert p.accepted
        assert np.allclose(p.delta, 0.0)

    def test_negative_dt_raises(self):
        with pytest.raises(ValueError):
            propose_step(cooling_series(), dt_years=-1.0, max_abs_step=50.0)

    def test_gate_config_passthrough(self):
        # loosening the curvature threshold flips a borderline refusal
        t = np.linspace(0, 10, 24)
        s = TrendSeries(times=t, values=260.0 - 1.0 * t - 0.03 * t**2)
        strict = propose_step(s, 10.0, 50.0, GateConfig(max_curvature_ratio=0.05))
        loose = propose_step(s, 10.0, 50.0, GateConfig(max_curvature_ratio=0.5))
        assert not strict.accepted
        assert loose.accepted
