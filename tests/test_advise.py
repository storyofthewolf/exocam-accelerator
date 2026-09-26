import numpy as np
import pytest

from synth import stefan_columns, truncate

from exocam_accelerate.advise import AdvisorConfig, advise, detect_jumps


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


def test_model_years_from_start_year():
    adv = advise(stefan_columns(years=50), "synth", start_year=101)
    assert adv.model_year == 150


class TestAfterAJump:
    def test_detects_jump_not_spinup(self):
        years = list(range(1, 101))
        cols = stefan_columns(h0=1.0, jump_year=61, factor=1.5)
        from exocam_accelerate.advise import _annual
        _, h = _annual(cols, "hi", "native")
        # early growth is >8 %/yr for years but is not a jump
        assert detect_jumps(years, h) == [61]

    def test_refuses_until_enough_post_jump_years(self):
        cols = stefan_columns(years=70, jump_year=61, factor=1.5)
        adv = advise(cols, "synth")
        assert adv.since_year == 61
        assert not adv.jump
        assert any("settled years since the jump" in r for r in adv.reasons)
        assert any("detected a jump" in w for w in adv.warnings)

    def test_fits_post_jump_native_data(self):
        cols = stefan_columns(years=61 + 2 + 20, jump_year=61, factor=1.5)
        adv = advise(cols, "synth")
        assert adv.jump, adv.reasons
        assert adv.which_used == "native"
        assert adv.window_years == 21      # every settled post-jump year
        a, b = adv.ice.fits["hyperbolic"].params
        assert b == pytest.approx(-70.0, rel=1e-2)

    def test_explicit_since(self):
        cols = stefan_columns(years=90)
        adv = advise(cols, "synth", AdvisorConfig(since_year=60))
        assert adv.since_year == 60
        assert adv.which_used == "native"
        assert adv.jump, adv.reasons
