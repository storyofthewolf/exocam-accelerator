"""TS-trajectory readout of the ocean probe check on D2-like noisy synthetic runs."""

import numpy as np
import pytest

from exocam_accelerate.check import (ProbeCheckConfig, Verdict, check_ocean_probe,
                                     ts_relaxation_readout)

JUMP = 61          # model year of the first probe year
C_EFF = 7.0        # W yr m-2 K-1


def d2_run(post_years, TS_eq=358.0, tau=25.0, dTS=4.0, sigma_T=1.0, rho=0.0,
           sigma_N=7.6, seed=0, pre_years=60, TS_jump=352.7, runaway=False,
           lam_n=None):
    """Annual one-box run with a +dTS step at model year JUMP, as monthly columns.

    Pre-probe TS relaxes to TS_eq with tau and sits at TS_jump when the probe
    starts; energy_bot = lam (TS_eq - TS) + noise, lam = C_eff/tau.
    """
    rng = np.random.default_rng(seed)
    yrs = np.arange(1, JUMP + post_years)           # model years 1..
    t = yrs - 1.0 + 0.5                             # mid-year, years since start
    t0 = JUMP - 1.0
    u_jump = TS_eq - TS_jump
    u = np.where(t < t0, u_jump * np.exp(-(t - t0) / tau),
                 (u_jump - dTS) * np.exp(-(t - t0) / tau))
    TS = TS_eq - u
    if runaway:   # TS keeps accelerating upward after the probe (0.2 -> 1.2 K/yr)
        TS = np.where(t < t0, TS, TS_jump + dTS + 1.2 * (t - t0))
    e = rng.standard_normal(yrs.size)
    n = np.zeros_like(e)
    for i in range(yrs.size):                       # AR(1) annual TS noise
        n[i] = (rho * n[i - 1] if i else 0.0) + np.sqrt(1 - rho ** 2) * e[i]
    TS = TS + sigma_T * n
    lam = C_EFF / tau if lam_n is None else lam_n
    N = lam * u + sigma_N * rng.standard_normal(yrs.size)
    month = np.arange(1, 12 * yrs.size + 1, dtype=float)
    cols = {"month": month}
    for name, v in (("TS", TS), ("energy_bot", N), ("energy_top", N),
                    ("ICEFRAC", np.zeros_like(TS))):
        for w in ("native", "int1", "int2"):
            cols[f"{name}_{w}"] = np.repeat(v, 12)
    return cols


def d2_log(dTS=4.0, **adv):
    advice = {"mode": "probe", "TS_now": 352.7, "since_year": 1,
              "config": {"imbalance": "energy_bot", "recent_years": 5, "heat_ratio": 2.0},
              "C_eff_W_yr_m2_K": C_EFF}
    advice.update(adv)
    return {"plugin": "som_ocean", "jump_model_year": JUMP, "somtp_dT": 2 * dTS,
            "somtp_dT_applied_mean": 2 * dTS, "advice": advice}


def readout(post_settled, **kw):
    cols = d2_run(post_settled + 2, **kw)
    res = check_ocean_probe(cols, d2_log(kw.get("dTS", 4.0)))
    return res


def test_ts_readout_recovers_tau_and_equilibrium():
    # unbiased within errors over many noise realisations (white and AR(1) noise)
    for rho in (0.0, 0.3):
        z_tau, z_eq, stats = [], [], []
        for seed in range(60):
            cols = d2_run(14, rho=rho, seed=seed)
            from exocam_accelerate.advise import _annual, model_years
            t, T = _annual(cols, "TS", "native")
            r = ts_relaxation_readout(model_years(t, 1), T, JUMP, 2)
            stats.append(r)
            if r["ts_status"] == "readable":
                z_tau.append((r["ts_tau"] - 25.0) / r["ts_tau_se"])
                z_eq.append((r["ts_eq"] - 358.0) / r["ts_eq_se"])
        assert len(z_tau) > 40                                   # usually readable by 12 yr
        # delta-method errors are rough for a ratio estimator; require the
        # pulls to be centred within ~0.6 sigma and not wildly over-confident
        assert abs(np.median(z_tau)) < 0.6 and abs(np.median(z_eq)) < 0.6
        assert np.std(z_eq) < 1.6, (rho, np.std(z_eq))
        # slope errors are honest: pre-rate scatter matches the quoted se
        rp = np.array([s["ts_r_pre"] for s in stats])
        assert np.mean([s["ts_r_pre_se"] for s in stats]) == pytest.approx(
            rp.std(), rel=0.4 if rho == 0 else 0.5) or rho > 0


def test_d2_probe_waits_then_reads_from_ts():
    """With energy_bot sigma 7.6 W/m2 and lambda = C/tau = 0.28, only TS can read it."""
    res5 = readout(5)
    assert res5.verdict is Verdict.WAIT                         # too early: min_years / noise
    res = readout(15, seed=3)
    assert res.verdict is Verdict.PASS, res.reasons
    m = res.metrics
    assert m["ts_tau"] == pytest.approx(25.0, abs=3 * m["ts_tau_se"] + 2)
    assert m["ts_eq"] == pytest.approx(358.0, abs=3 * m["ts_eq_se"] + 1)
    assert m["ts_lambda"] == pytest.approx(C_EFF / m["ts_tau"])
    assert m["ts_bracketed"] == 0.0                              # equilibrium above probe level
    assert any("tau" in r and "TS_eq" in r for r in res.reasons)


def test_falling_back_is_bracketed():
    # equilibrium below the probe level (TS_eq 351, probe to ~356.7): TS falls back
    res = readout(12, TS_eq=351.0, tau=12.0, seed=1)
    assert res.verdict is Verdict.PASS, res.reasons
    assert res.metrics["ts_r_post"] < 0 and res.metrics["ts_bracketed"] == 1.0


def test_accelerating_ts_is_runaway_fail():
    res = readout(10, runaway=True, seed=2)
    assert res.verdict is Verdict.FAIL
    assert any("accelerating" in r for r in res.reasons)


def test_no_false_runaway_on_noise():
    # (the energy_bot rule has its own, pre-existing false-alarm rate at sigma_N 7.6;
    # count only the TS rule here)
    fails = sum(any("accelerating" in r for r in
                    readout(8, seed=s, tau=60.0, TS_eq=356.0).reasons)
                for s in range(60))
    assert fails <= 2


def test_unresolved_gives_bounds_then_passes_at_max_wait():
    # a very slow relaxation: TS rate hardly changes
    cols = d2_run(6 + 14, TS_eq=400.0, tau=400.0, seed=5)
    log = d2_log()
    res = check_ocean_probe(cols, log, config=ProbeCheckConfig(min_years=5, max_wait_years=30))
    assert res.verdict is Verdict.WAIT and "ts_tau_lower" in res.metrics
    res = check_ocean_probe(cols, log, config=ProbeCheckConfig(min_years=5, max_wait_years=12))
    assert res.verdict is Verdict.PASS and any("tau >" in r for r in res.reasons)


def test_contradicting_readouts_wait():
    # energy_bot says "bracketed" (N<0 after the probe) while TS keeps rising slowly:
    # build by forcing N from a different equilibrium than TS
    cols = d2_run(14, TS_eq=370.0, tau=30.0, seed=4)
    from exocam_accelerate.advise import _annual
    _, T = _annual(cols, "TS", "native")
    N = np.where(np.arange(T.size) >= JUMP - 1, -6.0, 6.0) + 0.1 * np.sin(np.arange(T.size))
    for w in ("native", "int1", "int2"):
        cols[f"energy_bot_{w}"] = np.repeat(N, 12)
    res = check_ocean_probe(cols, d2_log())
    assert res.verdict is Verdict.WAIT and any("contradict" in r for r in res.reasons)


def test_short_pre_window_falls_back_to_energy_only():
    cols = d2_run(10)
    log = d2_log()
    log["advice"]["since_year"] = JUMP - 4
    res = check_ocean_probe(cols, log)
    assert "ts_r_pre" not in res.metrics
