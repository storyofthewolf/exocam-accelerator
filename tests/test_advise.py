import json

import numpy as np
import pytest

from synth import stefan_columns, truncate

from exocam_accelerate.advise import ADVICE_SCHEMA_VERSION, AdvisorConfig, advise, detect_jumps
from exocam_accelerate.cli import main
from exocam_accelerate.phase_space import PhaseGateConfig

_CAM_VARS = ("energy_top", "TS", "Tsfc", "ICEFRAC")
_CICE_VARS = ("hi", "qi")


def write_trend_files(directory, case_id, cols):
    """Serialize synth.stefan_columns() output as exocam-trend *_cam/_cice.txt."""
    month = cols["month"]
    n = month.size
    header_cols = {"cam": ["month"], "cice": ["month"]}
    rows = {"cam": [month], "cice": [month]}
    for var in _CAM_VARS + _CICE_VARS:
        comp = "cam" if var in _CAM_VARS else "cice"
        for which in ("native", "int1", "int2"):
            key = f"{var}_{which}" if which != "int1" else f"{var}_native"
            v = cols.get(f"{var}_{which}", cols[f"{var}_native"])
            header_cols[comp].append(f"{var}_{which}")
            rows[comp].append(v)
    span = f"0001-01-{1 + (n - 1) // 12:04d}-12"
    for comp in ("cam", "cice"):
        path = directory / f"{case_id}_{span}_{comp}.txt"
        data = np.column_stack(rows[comp])
        with path.open("w") as f:
            f.write("  ".join(header_cols[comp]) + "\n")
            np.savetxt(f, data)
    return directory


def test_recommends_factor_from_conduction_law():
    cols = stefan_columns()
    adv = advise(cols, "synth", AdvisorConfig(max_ice_factor=2.0))
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


def test_missing_icefrac_refuses():
    cols = {k: v for k, v in stefan_columns().items() if not k.startswith("ICEFRAC")}
    adv = advise(cols, "synth")
    assert not adv.jump
    assert any("ICEFRAC" in r for r in adv.reasons)


def test_missing_qi_refuses():
    cols = {k: v for k, v in stefan_columns().items() if not k.startswith("qi_")}
    adv = advise(cols, "synth")
    assert not adv.jump
    assert any("qi" in r for r in adv.reasons)


def test_qi_hi_drift_beyond_threshold_refuses():
    cols = stefan_columns()
    t = np.arange(1, cols["month"].size + 1) / 12.0
    # qi/hi ratio should be constant by default (qi ~ hi); force it to drift.
    cols["qi_native"] = cols["qi_native"] * (1.0 + 0.5 * t / t[-1])
    cols["qi_int2"] = cols["qi_int2"] * (1.0 + 0.5 * t / t[-1])
    adv = advise(cols, "synth")
    assert not adv.jump
    assert any("qi/hi drifted" in r for r in adv.reasons)


def test_rejected_temperature_reference_refuses():
    cols = stefan_columns()
    # Break TS's linear relation to N without touching the ice law, so the
    # ice extrapolation itself would otherwise accept.
    rng = np.random.default_rng(0)
    noise = 50.0 * rng.standard_normal(cols["TS_native"].size)
    cols["TS_native"] = cols["TS_native"] + noise
    cols["TS_int2"] = cols["TS_int2"] + noise
    adv = advise(cols, "synth")
    assert not adv.jump
    assert any("TS(N) reference rejected" in r for r in adv.reasons)


def test_explicit_n_target_must_advance_without_reversing():
    # a=0.1, N_now_fit ~ -4.57 for the default stefan_columns(); a target more
    # negative than N_now_fit moves further from equilibrium, not toward it.
    adv = advise(stefan_columns(), "synth", AdvisorConfig(N_target=-10.0))
    assert not adv.jump
    assert any("does not lie strictly between" in r for r in adv.reasons)


def test_to_dict_has_schema_version_and_provenance():
    adv = advise(stefan_columns(), "synth", provenance={"a.txt": "deadbeef"})
    d = adv.to_dict()
    assert d["schema_version"] == ADVICE_SCHEMA_VERSION
    assert d["provenance"] == {"a.txt": "deadbeef"}


class TestConfigDomainValidation:
    @pytest.mark.parametrize("kwargs", [
        dict(n_fraction=0.0), dict(n_fraction=1.0), dict(n_fraction=-0.1),
        dict(max_ice_factor=0.5), dict(max_ice_factor=2.1),
        dict(window_years=0.0), dict(window_years=-5.0),
        dict(settle_years=0), dict(post_jump_min_years=0),
    ])
    def test_invalid_advisor_config_rejected(self, kwargs):
        with pytest.raises(ValueError):
            AdvisorConfig(**kwargs)

    def test_invalid_gate_config_rejected(self):
        with pytest.raises(ValueError):
            PhaseGateConfig(max_extrapolation_ratio=0.0)
        with pytest.raises(ValueError):
            PhaseGateConfig(max_extrapolation_ratio=-1.0)


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


class TestAdviseCLI:
    def test_json_output_carries_schema_version_and_provenance(self, tmp_path):
        write_trend_files(tmp_path, "synth", stefan_columns())
        out = tmp_path / "advice.json"
        assert main(["advise", str(tmp_path), "synth", "--max-ice-factor", "2.0",
                     "--json", str(out)]) == 0
        adv = json.loads(out.read_text())
        assert adv["schema_version"] == ADVICE_SCHEMA_VERSION
        assert set(adv["provenance"]) == {"synth_0001-01-0100-12_cam.txt",
                                          "synth_0001-01-0100-12_cice.txt"}

    def test_invalid_parameter_domain_rejected(self, tmp_path):
        write_trend_files(tmp_path, "synth", stefan_columns())
        assert main(["advise", str(tmp_path), "synth", "--n-fraction", "0.0"]) == 2
