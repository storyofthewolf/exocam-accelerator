import numpy as np
import pytest

from exocam_accelerate.plugins import PLUGIN_REGISTRY
from exocam_accelerate.som_ocean import (MAX_DT, TK_FRZ_SW, PatternConfig,
                                         SomOceanPlugin, area_mean, check_dT,
                                         pattern_weights, to_grid)


def test_registered():
    assert "som_ocean" in PLUGIN_REGISTRY


def test_to_grid_is_lon_fastest():
    v = np.arange(12)
    g = to_grid(v, 3, 4)
    assert g[1, 0] == 4 and g[0, 3] == 3
    with pytest.raises(ValueError):
        to_grid(v, 4, 4)


def test_uniform_shift_and_constraints():
    p = SomOceanPlugin()
    before = {"somtp": np.array([300.0, 310.0, TK_FRZ_SW, 290.0])}
    after = p.apply_delta(before, 2.0)
    out, rep = p.enforce_constraints(before, after,
                                     ocean=np.array([True, True, True, False]))
    np.testing.assert_allclose(out["somtp"], [302.0, 312.0, TK_FRZ_SW, 290.0])
    msg = rep.adjustments["somtp"]
    assert "near-freezing" in msg and "non-ocean" in msg
    assert before["somtp"][0] == 300.0          # inputs untouched


def test_cooling_is_clamped_at_freezing():
    p = SomOceanPlugin()
    before = {"somtp": np.array([272.0, 300.0])}
    out, rep = p.enforce_constraints(before, p.apply_delta(before, -5.0))
    assert out["somtp"][0] == pytest.approx(TK_FRZ_SW)
    assert "clamped" in rep.adjustments["somtp"]


def test_bounds_and_sanity():
    with pytest.raises(ValueError):
        check_dT(MAX_DT + 1)
    with pytest.raises(ValueError):
        check_dT(np.array([1.0, np.nan]))
    p = SomOceanPlugin()
    with pytest.raises(ValueError):          # Celsius, not Kelvin
        p.enforce_constraints({"somtp": np.array([25.0])}, {"somtp": np.array([26.0])})
    with pytest.raises(ValueError):
        p.enforce_constraints({"somtp": np.array([np.nan])}, {"somtp": np.array([1.0])})
    with pytest.raises(ValueError):
        p.apply_delta({"somtp": np.zeros(4) + 300}, np.ones(3))


def test_pattern_weights_area_mean_one_and_bounds():
    nj, ni = 6, 8
    lat = np.linspace(-75, 75, nj)[:, None] * np.ones((1, ni))
    area = np.cos(np.radians(lat))
    old = 330.0 + np.zeros((nj, ni))
    rate = 0.3 + 0.2 * np.cos(np.radians(lat))          # tropics warm faster
    now = old + 10.0 * rate
    pat = pattern_weights(old, now, 10.0, area, config=PatternConfig(smooth=0))
    assert area_mean(pat.weight, area) == pytest.approx(1.0)
    np.testing.assert_allclose(pat.weight, rate / area_mean(rate, area), rtol=1e-9)
    assert pat.mean_rate == pytest.approx(area_mean(rate, area))
    s = pat.summary()
    assert s["weight_min"] < 1 < s["weight_max"]


def test_pattern_clips_and_renormalizes():
    nj, ni = 4, 4
    area = np.ones((nj, ni))
    old = np.full((nj, ni), 330.0)
    now = old + 1.0
    now[0, 0] += 30.0                                     # one hot outlier
    pat = pattern_weights(old, now, 10.0, area,
                          config=PatternConfig(smooth=0, max_weight=2.0))
    assert pat.weight.max() <= 2.0 + 1e-9
    assert area_mean(pat.weight, area) == pytest.approx(1.0, abs=1e-6)


def test_pattern_excludes_frozen_and_masked_cells():
    area = np.ones((3, 3))
    old = np.full((3, 3), 300.0)
    old[0, 0] = TK_FRZ_SW
    now = old + 2.0
    now[0, 0] = TK_FRZ_SW
    ocean = np.ones((3, 3), bool)
    ocean[2, 2] = False
    pat = pattern_weights(old, now, 10.0, area, ocean, PatternConfig(smooth=1))
    assert pat.weight[0, 0] == 0 and pat.weight[2, 2] == 0
    assert area_mean(pat.weight, pat.area) == pytest.approx(1.0)


def test_pattern_refuses_noise():
    area = np.ones((3, 3))
    rng = np.random.default_rng(0)
    old = 330 + rng.normal(0, 0.1, (3, 3))
    with pytest.raises(ValueError, match="noise"):
        pattern_weights(old, old + 0.01, 10.0, area)
