import numpy as np
import pytest

from exocam_accelerate.aqua_ice import AquaIcePlugin, check_factor
from exocam_accelerate.plugins import PLUGIN_REGISTRY


def fields():
    rng = np.random.default_rng(1)
    vicen = rng.uniform(0, 30, (5, 4, 6))
    vicen[:, 0, :] = 0.0                       # ice-free row
    eicen = -rng.uniform(1e8, 1e9, (20, 4, 6))
    eicen[:, 0, :] = 0.0
    vsnon = rng.uniform(0, 0.2, (5, 4, 6))
    esnon = -rng.uniform(1e6, 1e7, (5, 4, 6))
    return dict(vicen=vicen, eicen=eicen, vsnon=vsnon, esnon=esnon)


def test_registered():
    assert isinstance(PLUGIN_REGISTRY["aqua_ice"], AquaIcePlugin)


def test_scales_ice_only_by_default():
    f = fields()
    out = AquaIcePlugin().apply_delta(f, 1.4)
    assert set(out) == {"vicen", "eicen"}
    np.testing.assert_allclose(out["vicen"], 1.4 * f["vicen"])
    np.testing.assert_allclose(out["eicen"], 1.4 * f["eicen"])
    # enthalpy per unit volume unchanged (per category sums)
    assert f["vicen"] is not out["vicen"]


def test_separate_snow_factor():
    f = fields()
    out = AquaIcePlugin().apply_delta(f, (1.4, 1.1))
    np.testing.assert_allclose(out["vsnon"], 1.1 * f["vsnon"])
    np.testing.assert_allclose(out["esnon"], 1.1 * f["esnon"])


def test_constraints_clean_state_untouched():
    f = fields()
    p = AquaIcePlugin()
    after = p.apply_delta(f, 1.3)
    fixed, report = p.enforce_constraints(f, after)
    np.testing.assert_array_equal(fixed["vicen"], after["vicen"])
    assert "x1.3000" in report.adjustments["vicen"]
    assert "clamped" not in report.adjustments["vicen"]


def test_constraints_clamp_bad_input():
    f = fields()
    f["vicen"][1, 1, 1] = -1.0
    f["eicen"][2, 2, 2] = 5.0
    p = AquaIcePlugin()
    fixed, report = p.enforce_constraints(f, p.apply_delta(f, 1.2))
    assert fixed["vicen"][1, 1, 1] == 0.0
    assert fixed["eicen"][2, 2, 2] == 0.0
    assert "clamped 1" in report.adjustments["vicen"]
    assert "clamped 1" in report.adjustments["eicen"]


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), 2.5, 0.3])
def test_factor_bounds(bad):
    with pytest.raises(ValueError):
        check_factor(bad)
