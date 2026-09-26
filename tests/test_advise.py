import numpy as np
import pytest

from exocam_accelerate.advise import AdvisorConfig, advise


def stefan_columns(years=100, a=0.1, b=-70.0, h0=15.0, k=6.0,
                   icefrac=0.8, icefrac_drift=0.0):
    """Synthetic trend columns obeying conduction-limited ice growth.

    hi(t) = sqrt(h0^2 + 2 k t),  N = a + b/hi,  TS linear in N,  qi ∝ hi.
    int2 is set equal to native (smoothing is not under test here).
    """
    month = np.arange(1, 12 * years + 1, dtype=float)
    t = month / 12.0
    hi = np.sqrt(h0**2 + 2 * k * t)
    N = a + b / hi
    series = {
        "hi": hi,
        "energy_top": N,
        "TS": 210.0 + 1.5 * N,
        "Tsfc": -60.0 + 1.5 * N,
        "qi": -3.0e20 * hi,
        "ICEFRAC": icefrac + icefrac_drift * t / years,
    }
    cols = {"month": month}
    for name, v in series.items():
        cols[f"{name}_native"] = v
        cols[f"{name}_int2"] = v
    return cols


def test_recommends_factor_from_conduction_law():
    cols = stefan_columns()
    adv = advise(cols, "synth", AdvisorConfig(max_ice_factor=10.0))
    assert adv.jump, adv.reasons
    a, b = adv.ice.fits["hyperbolic"].params
    # annual averaging of a nonlinear law biases the fit very slightly
    assert a == pytest.approx(0.1, abs=1e-3)
    assert b == pytest.approx(-70.0, rel=1e-3)
    # target removes half the (fitted) imbalance
    assert adv.N_target == pytest.approx(0.5 * adv.N_now_fit)
    h_now = adv.ice.X_now
    expected = (b / (adv.N_target - a)) / h_now
    assert adv.ice_factor == pytest.approx(expected, rel=1e-6)
    assert not adv.ice_factor_clipped
    assert adv.N_after == pytest.approx(adv.N_target, rel=1e-6)
    # TS reference follows the linear law at N_after
    assert adv.temperatures["TS"].prediction == pytest.approx(
        210.0 + 1.5 * adv.N_after, abs=1e-6)
    # Stefan: years skipped from the h^2 slope (2k per year)
    h_new = h_now * adv.ice_factor
    assert adv.years_skipped == pytest.approx((h_new**2 - h_now**2) / 12.0, rel=1e-3)


def test_factor_is_clipped():
    adv = advise(stefan_columns(), "synth", AdvisorConfig(max_ice_factor=1.2))
    assert adv.jump
    assert adv.ice_factor == pytest.approx(1.2)
    assert adv.ice_factor_clipped
    assert adv.ice_factor_raw > 1.2
    # N_after reflects the clipped jump, not the unreachable target
    assert adv.N_after < adv.N_target


def test_auto_window_prefers_longest_settled():
    adv = advise(stefan_columns(), "synth")
    assert adv.window_years == 40.0


def test_refuses_when_ice_edge_moving():
    cols = stefan_columns(icefrac_drift=0.5)
    adv = advise(cols, "synth")
    assert not adv.jump
    assert any("ice edge" in r for r in adv.reasons)


def test_refuses_target_past_asymptote():
    adv = advise(stefan_columns(), "synth", AdvisorConfig(N_target=0.5))
    assert not adv.jump


def test_missing_hi():
    cols = {k: v for k, v in stefan_columns().items() if not k.startswith("hi_")}
    adv = advise(cols, "synth")
    assert not adv.jump
    assert adv.ice is None


def test_to_dict_is_json_serialisable():
    import json
    adv = advise(stefan_columns(), "synth")
    d = json.loads(json.dumps(adv.to_dict()))
    assert d["ice_factor"] == pytest.approx(adv.ice_factor)
    assert "TS" in d["expected_after_jump"]
