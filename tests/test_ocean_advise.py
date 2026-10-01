import numpy as np
import pytest
from synth import gregory_columns, truncate

from exocam_accelerate.ocean_advise import (OCEAN_ADVICE_SCHEMA_VERSION,
                                            OceanAdvisorConfig, advise_ocean,
                                            detect_ts_jumps, pattern_advice)

C0, C1, TAU = 350.0, -1.2, 12.0


def test_recovers_gregory_line_and_target():
    a = advise_ocean(truncate(gregory_columns(), 30), "synth")
    assert a.jump, a.reasons
    c0, c1 = a.gregory.fits["linear"].params
    assert c0 == pytest.approx(C0, abs=1e-6)
    assert c1 == pytest.approx(C1, abs=1e-6)
    # half the imbalance removed: TS goes halfway to c0
    assert a.TS_after == pytest.approx(a.TS_now + 0.5 * (C0 - a.TS_now), abs=1e-6)
    assert a.N_after == pytest.approx(0.5 * a.N_now_fit, rel=1e-6)
    assert a.somtp_dT == pytest.approx(a.TS_after - a.TS_now)


def test_heat_capacity_and_years_skipped_are_one_box_exact():
    a = advise_ocean(truncate(gregory_columns(), 30), "synth")
    # C = -tau/c1; skipping half the imbalance takes tau*ln 2
    assert a.C_eff == pytest.approx(-TAU / C1, rel=0.02)
    assert a.tau_years == pytest.approx(TAU, rel=0.02)
    assert a.years_skipped == pytest.approx(TAU * np.log(2), rel=0.02)


def test_default_coordinate_is_energy_bot_and_gap_reported():
    # atmospheric storage: energy_top = energy_bot + C_atm dTS/dt
    cols = truncate(gregory_columns(c_atm=20.0), 30)
    a = advise_ocean(cols, "synth")
    assert a.config.imbalance == "energy_bot"
    assert a.gap["mean"] > 0
    # the equilibrium from the surface line is exact; the TOA line is bent
    assert a.gregory.fits["linear"].params[0] == pytest.approx(C0, abs=1e-6)
    top = advise_ocean(cols, "synth", OceanAdvisorConfig(imbalance="energy_top"))
    assert any("storage" in w for w in top.warnings) or not top.jump


def test_heat_ratio_measured_from_storage():
    # C_ocean = -tau/c1 = 10; the atmosphere stores 20 more per K -> ratio 3
    a = advise_ocean(truncate(gregory_columns(c_atm=20.0), 30), "synth")
    assert a.heat["measured"]
    assert a.heat["C_ocean"] == pytest.approx(10.0, rel=0.03)
    assert a.heat["ratio"] == pytest.approx(3.0, rel=0.05)
    assert a.somtp_dT == pytest.approx(a.heat["used"] * (a.TS_after - a.TS_now))
    assert a.to_dict()["config"]["heat_ratio"] == pytest.approx(a.heat["used"])
    big = advise_ocean(truncate(gregory_columns(c_atm=60.0), 30), "synth")
    assert big.heat["used"] == pytest.approx(4.0)         # clipped


def test_heat_ratio_scales_increment_not_landing():
    cols = truncate(gregory_columns(), 30)
    a = advise_ocean(cols, "synth")
    assert a.heat["used"] == pytest.approx(1.0, abs=0.02)
    b = advise_ocean(cols, "synth", OceanAdvisorConfig(heat_ratio=1.5))
    assert b.somtp_dT == pytest.approx(1.5 * a.somtp_dT, rel=0.02)
    assert b.TS_after == pytest.approx(a.TS_after)


def test_clip():
    a = advise_ocean(truncate(gregory_columns(T0=250.0), 30), "synth",
                     OceanAdvisorConfig(max_dT=2.0))
    assert a.somtp_dT == pytest.approx(2.0)
    assert a.clipped
    assert a.TS_after == pytest.approx(a.TS_now + 2.0)


def test_cooling_run_gets_negative_increment():
    a = advise_ocean(truncate(gregory_columns(T0=400.0), 30), "synth")
    assert a.jump and a.somtp_dT < 0


def test_refuses_with_ice():
    a = advise_ocean(truncate(gregory_columns(icefrac=0.05), 30), "synth")
    assert not a.jump
    assert any("ice-free" in r for r in a.reasons)


def test_early_icy_years_are_skipped_not_refused():
    cols = truncate(gregory_columns(), 30)
    cols["ICEFRAC_native"] = np.where(cols["month"] <= 12 * 8, 0.05, 0.0)
    a = advise_ocean(cols, "synth")
    assert a.jump, a.reasons


def test_refuses_runaway_like_relation():
    cols = truncate(gregory_columns(), 30)
    for w in ("native", "int1", "int2"):   # N grows with TS: no equilibrium
        cols[f"energy_bot_{w}"] = -cols[f"energy_bot_{w}"]
    a = advise_ocean(cols, "synth")
    assert not a.jump


def test_refuses_when_current_state_is_off_the_line():
    # the last years sit well off the relation the window shows (steepening)
    cols = truncate(gregory_columns(), 40)
    late = cols["month"] > 12 * 35
    for w in ("native", "int1", "int2"):
        cols[f"energy_bot_{w}"] = np.where(late, cols[f"energy_bot_{w}"] + 3.0,
                                           cols[f"energy_bot_{w}"])
    a = advise_ocean(cols, "synth", OceanAdvisorConfig(window_years=30))
    assert not a.jump
    assert any("off the line" in r for r in a.reasons)


def test_near_equilibrium_not_worth_it():
    a = advise_ocean(gregory_columns(years=120, noise=0.3), "synth")
    assert not a.jump
    assert any("not worth" in r for r in a.reasons)


def test_nan_in_latest_year_refuses():
    cols = truncate(gregory_columns(), 30)
    cols["TS_native"] = cols["TS_native"].copy()
    cols["TS_native"][-12:] = np.nan
    a = advise_ocean(cols, "synth")
    assert not a.jump


def test_post_jump_fit_uses_native_after_since():
    cols = gregory_columns(years=60, tau=30.0, jump_year=31, dT=5.0)
    a = advise_ocean(cols, "synth", OceanAdvisorConfig(since_year=31))
    assert a.since_year == 31 and a.which_used == "native"
    assert a.jump, a.reasons
    assert a.gregory.fits["linear"].params[0] == pytest.approx(C0, abs=1e-6)
    few = advise_ocean(truncate(cols, 40), "synth", OceanAdvisorConfig(since_year=31))
    assert not few.jump and any("settled years" in r for r in few.reasons)


def test_detector_finds_a_clean_step_but_is_off_by_default():
    cols = gregory_columns(years=60, jump_year=41, dT=4.0, noise=0.05, seed=3)
    from exocam_accelerate.advise import _annual, model_years
    t, T = _annual(cols, "TS", "native")
    assert 41 in detect_ts_jumps(model_years(t, 1), T)
    assert advise_ocean(cols, "synth").since_year is None


def test_serialization_and_pattern_advice():
    d = advise_ocean(truncate(gregory_columns(), 30), "synth").to_dict()
    assert d["schema_version"] == OCEAN_ADVICE_SCHEMA_VERSION
    assert d["plugin"] == "som_ocean"
    assert d["expected_after_jump"]["TS"] == d["TS_after"]
    p = pattern_advice(d, {"weight_min": 0.5})
    assert p["schema_version"] == "ocean-2" and p["somtp_dT"] == d["somtp_dT"]
    with pytest.raises(ValueError):
        pattern_advice(p, {})


def test_config_validation():
    with pytest.raises(ValueError):
        OceanAdvisorConfig(max_dT=40.0)
    with pytest.raises(ValueError):
        OceanAdvisorConfig(heat_ratio=0.5)
    with pytest.raises(ValueError):
        OceanAdvisorConfig(imbalance="FLNT")


# ---------------------------------------------------------------- probe mode
from exocam_accelerate.ocean_advise import OCEAN_PROBE_SCHEMA_VERSION, ProbeConfig, probe_ocean  # noqa: E402


def test_probe_steps_by_recent_trend():
    cols = truncate(gregory_columns(tau=40.0, noise=0.02), 40)
    d = probe_ocean(cols, "synth")
    assert d["schema_version"] == OCEAN_PROBE_SCHEMA_VERSION and d["mode"] == "probe"
    assert d["dTS"] == pytest.approx(d["trend_K_per_yr"] * 15.0, rel=0.05)
    s = probe_ocean(truncate(gregory_columns(tau=40.0, noise=0.02, c_atm=20.0), 40), "synth")
    assert s["heat"]["used"] == pytest.approx(1 + 20.0 / (40.0 / 1.2), rel=0.05)
    assert s["somtp_dT"] == pytest.approx(
        s["heat"]["used"] * s["dTS"])
    assert d["TS_after"] == pytest.approx(d["TS_now"] + d["somtp_dT"])
    assert d["lambda_detectable"] > 0


def test_probe_clip_explicit_and_refusals():
    cols = truncate(gregory_columns(tau=40.0, noise=0.02), 40)
    assert probe_ocean(cols, "s", ProbeConfig(max_dT=1.0))["somtp_dT"] == pytest.approx(1.0)
    d = probe_ocean(cols, "s", ProbeConfig(probe_dT=-2.0))
    assert d["dTS"] == pytest.approx(-2.0) and d["sizing"] == "explicit"
    assert probe_ocean(truncate(gregory_columns(icefrac=0.1), 40), "s")["somtp_dT"] is None
    # a flat, noisy run: the probe would drown in the noise
    flat = gregory_columns(years=150, noise=0.5)
    d = probe_ocean(flat, "s", ProbeConfig(probe_dT=0.3))
    assert d["somtp_dT"] is None and any("sigmas" in r for r in d["reasons"])
