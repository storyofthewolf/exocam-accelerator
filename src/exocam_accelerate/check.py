"""Post-jump consistency check: did the run continue along the predicted path?

Pure computation. Implements the design of ``consistency.py`` for the aqua_ice
jump. Inputs are the case's exocam-trend columns after the jump and the jump
log (the ``<cice.r>.accel.json`` sidecar ``restart.apply_ice_jump`` writes),
which holds the advice the jump was sized from.

The reference is the advice's own fitted relations, evaluated at the ice the
run actually has each year — not a single number:

* landed:   first post-jump annual hi ~ hi_before x factor (else the edit
            never reached the model: wrong file, stale rpointer, rollback);
* on-law:   N_obs ~ a + b/hi_obs   (the conduction law the jump relied on;
            plus ``law_offset_after`` for a tapered jump);
* TS:       TS_obs ~ c0 + c1*N_law (the Gregory relation it relied on);
  both measured relative to the run's own offset from those relations over the
  last pre-jump years, so slow drift of the fit is not blamed on the jump;
* holding:  hi is not melting back.

Verdicts: PASS, WAIT (not enough settled years yet, nothing alarming), FAIL
(roll back). The tolerances come from the end-to-end hindcast (p90 |N_after
error| 0.34 W/m2, |TS| 0.44 K) widened by the interannual noise of an n-year
native mean.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .advise import ICE_VOLUME_VAR, _annual, _have, model_years


class Verdict(enum.Enum):
    PASS = "PASS"
    WAIT = "WAIT"
    FAIL = "FAIL"

    @property
    def exit_code(self) -> int:
        return {"PASS": 0, "WAIT": 10, "FAIL": 20}[self.value]


@dataclass(frozen=True)
class CheckConfig:
    settle_years: int = 2         # adjustment years excluded from the test
    min_years: int = 3            # settled years needed before PASS
    tol_N: float = 0.4            # W/m2, floor on the N tolerance
    tol_TS: float = 0.75          # K, floor on the TS tolerance (secondary to N)
    land_tol: float = 0.15        # |hi_first / hi_expected - 1|
    #: an early FAIL needs a deviation this many tolerances out
    early_fail_factor: float = 3.0
    #: pre-jump years whose offset from the fitted relations is subtracted
    baseline_years: int = 5

    def __post_init__(self) -> None:
        if self.settle_years <= 0:
            raise ValueError(f"settle_years must be positive, got {self.settle_years!r}")
        if self.min_years <= 0:
            raise ValueError(f"min_years must be positive, got {self.min_years!r}")
        if self.baseline_years <= 0:
            raise ValueError(f"baseline_years must be positive, got "
                             f"{self.baseline_years!r}")


@dataclass(frozen=True)
class CheckResult:
    verdict: Verdict
    reasons: List[str]
    jump_year: int
    years_after: int              # post-jump years available
    settled_years: int
    metrics: Dict[str, float] = field(default_factory=dict)
    #: flags reported alongside the verdict without deciding it (e.g.
    #: "runaway greenhouse suspected"): what to do about them is the user's call
    warnings: List[str] = field(default_factory=list)


def _noise(x) -> float:
    """Interannual std estimated from year-to-year differences (trend-robust)."""
    x = np.asarray(x, dtype=float)
    if x.size < 3:
        return 0.0
    return float(np.std(np.diff(x)) / np.sqrt(2.0))


def check_jump(columns: Dict[str, np.ndarray], jump_log: dict,
               start_year: int = 1,
               config: CheckConfig = CheckConfig()) -> CheckResult:
    """Score the run after a jump against the advice recorded in ``jump_log``."""
    advice = jump_log.get("advice") or {}
    jump_year = int(jump_log["jump_model_year"])
    # a tapered jump (schema 2) scales cells unequally: the global mean moves
    # by the effective factor, and the run should sit law_offset_after off the
    # global law (taper.py) — zero for a uniform jump
    factor = float(advice.get("effective_factor") or jump_log["ice_factor"])
    offset = float(advice.get("law_offset_after") or 0.0)
    ice = advice.get("ice") or {}
    try:
        a, b = ice["fits"]["hyperbolic"]["params"]
        h_before = float(ice["X_now"])
    except (KeyError, TypeError, ValueError):
        raise ValueError("jump log carries no conduction-law fit (advice.ice.fits."
                         "hyperbolic); was the jump made with --ice-factor alone?")
    imbalance = (advice.get("config") or {}).get("imbalance", "energy_top")

    t, N = _annual(columns, imbalance, "native")
    _, h = _annual(columns, ICE_VOLUME_VAR, "native")
    years = model_years(t, start_year)
    post = years >= jump_year
    settled = years >= jump_year + config.settle_years
    n_post, n_set = int(post.sum()), int(settled.sum())
    reasons: List[str] = []
    metrics: Dict[str, float] = {}

    def done(verdict):
        return CheckResult(verdict, reasons, jump_year, n_post, n_set, metrics)

    if n_post == 0:
        reasons.append(f"no complete model year since the jump at year {jump_year} "
                       f"in the trend data (latest is {int(years[-1])})")
        return done(Verdict.WAIT)

    # landed?
    h_expected = h_before * factor
    h_first = float(h[post][0])
    land_err = h_first / h_expected - 1.0
    metrics.update(hi_expected=h_expected, hi_first=h_first, land_error=land_err)
    if abs(land_err) > config.land_tol:
        reasons.append(f"jump did not land: first post-jump hi {h_first:.2f} m vs "
                       f"{h_expected:.2f} expected ({land_err:+.0%}); was the "
                       f"edited cice.r the one the run read?")
        return done(Verdict.FAIL)

    pre = (years < jump_year) & (years >= jump_year - 10)
    sig_N = _noise(N[pre]) if pre.sum() >= 3 else _noise(N[post])

    if n_set == 0:
        reasons.append(f"jump landed; in the {config.settle_years}-yr adjustment "
                       f"period ({n_post} yr so far)")
        return done(Verdict.WAIT)

    # on the conduction law? Measured relative to how far the run already sat
    # off the fitted law just before the jump (the law drifts slowly; that
    # offset is not the jump's doing).
    last = (years < jump_year) & (years >= jump_year - config.baseline_years)
    bias_N = float(np.mean(N[last] - (a + b / h[last]))) if last.any() else 0.0
    N_law = a + b / h[settled] + offset
    dN = float(np.mean(N[settled] - N_law)) - bias_N
    tol_N = max(config.tol_N, 2.0 * sig_N / np.sqrt(n_set))
    metrics.update(N_obs=float(N[settled].mean()), N_law=float(N_law.mean()),
                   dN=dN, tol_N=tol_N, N_bias_before=bias_N)

    # TS on the Gregory relation? Only when the advisor's TS(N) reference was
    # actually accepted (feasibility-review finding 2/1): a rejected or
    # missing reference is a defect in the advice itself, not something more
    # settled years can resolve, so it must not be used, and it forces the
    # verdict away from PASS below (`ts_ok`).
    dT = tol_T = None
    ts = (advice.get("temperatures") or {}).get("TS") or {}
    ts_ok = bool(ts.get("accepted"))
    lin = (ts.get("fits") or {}).get("linear")
    if ts_ok and lin and f"TS_native" in columns:
        c0, c1 = lin["params"]
        _, T = _annual(columns, "TS", "native")
        bias_T = (float(np.mean(T[last] - (c0 + c1 * (a + b / h[last]))))
                  if last.any() else 0.0)
        dT = float(np.mean(T[settled] - (c0 + c1 * N_law))) - bias_T
        sig_T = _noise(T[pre]) if pre.sum() >= 3 else _noise(T[post])
        tol_T = max(config.tol_TS, 2.0 * sig_T / np.sqrt(n_set))
        metrics.update(TS_obs=float(T[settled].mean()), dTS=dT, tol_TS=tol_T)

    # holding? (only meaningful with a few years)
    if n_post >= 3:
        slope = float(np.polyfit(t[post], h[post], 1)[0])
        metrics["hi_trend_m_per_yr"] = slope
        if slope < 0:
            reasons.append(f"ice is melting back after the jump "
                           f"({slope:+.2f} m/yr)")

    bad_N = abs(dN) > tol_N
    bad_T = dT is not None and abs(dT) > tol_T
    if bad_N:
        reasons.append(f"N is off the conduction law by {dN:+.2f} W/m2 "
                       f"(tolerance {tol_N:.2f})")
    if bad_T:
        reasons.append(f"TS is off the Gregory relation by {dT:+.2f} K "
                       f"(tolerance {tol_T:.2f})")
    if not ts_ok:
        reasons.append("no accepted TS(N) reference in the jump log: the standard "
                       "aqua_ice protocol cannot verify TS, so this jump cannot PASS")
    elif dT is None:
        reasons.append("no TS reference in the jump log: TS not checked")

    melting = "hi_trend_m_per_yr" in metrics and metrics["hi_trend_m_per_yr"] < 0
    if n_set < config.min_years:
        k = config.early_fail_factor
        if abs(dN) > k * tol_N or (dT is not None and abs(dT) > k * tol_T):
            return done(Verdict.FAIL)
        reasons.append(f"{n_set} of {config.min_years} settled years so far")
        return done(Verdict.WAIT)
    # A missing/rejected TS reference is a defect in the advice, not in the
    # run: once enough settled years exist to otherwise PASS, treat it as
    # FAIL (not WAIT) — waiting longer cannot fix a reference that was never
    # accepted, and the standard protocol requires it (feasibility-review
    # finding 2).
    if bad_N or bad_T or melting or not ts_ok:
        return done(Verdict.FAIL)
    return done(Verdict.PASS)


# ---------------------------------------------------------------------------
# som_ocean jump (docn.r somtp; ocean_advise.py)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OceanCheckConfig:
    settle_years: int = 2
    min_years: int = 3
    #: W/m2 floor on the tolerance for N off the Gregory line (energy_bot is
    #: noisier than energy_top: ~1-3 W/m2 interannual in the hot runs)
    tol_N: float = 0.75
    #: the jump "did not land" when the first post-jump year moved TS by less
    #: than this fraction of the expected change (after the pre-jump drift)...
    min_land_fraction: float = 0.3
    #: ...and the shortfall is larger than this many interannual TS sigmas
    land_sigmas: float = 2.0
    early_fail_factor: float = 3.0
    baseline_years: int = 5

    def __post_init__(self) -> None:
        for name in ("settle_years", "min_years", "baseline_years"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")


def check_ocean_jump(columns: Dict[str, np.ndarray], jump_log: dict,
                     start_year: int = 1,
                     config: OceanCheckConfig = OceanCheckConfig()) -> CheckResult:
    """Score the run after a somtp jump against the advice in ``jump_log``.

    * landed:  the first post-jump annual TS moved by the expected
               ``somtp_dT / heat_ratio`` (net of the pre-jump drift); the
               ratio observed is reported as ``heat_ratio_implied`` —
               the measurement the next jump's ``--heat-ratio`` needs;
    * on-line: settled-year N ~ (TS - c0)/c1, the Gregory line the jump relied
               on, relative to the run's offset from it just before the jump;
    * sanity:  non-finite data in the post-jump years is a FAIL, never a PASS.
    """
    advice = jump_log.get("advice") or {}
    if advice.get("mode") == "probe":
        return check_ocean_probe(columns, jump_log, start_year, ProbeCheckConfig(
            settle_years=config.settle_years))
    jump_year = int(jump_log["jump_model_year"])
    greg = advice.get("gregory") or {}
    try:
        c0, c1 = greg["fits"]["linear"]["params"]
        TS_before = float(advice["TS_now"])
    except (KeyError, TypeError, ValueError):
        raise ValueError("jump log carries no Gregory fit (advice.gregory.fits.linear); "
                         "was the jump made with --delta-t alone?")
    if not greg.get("accepted"):
        raise ValueError("the jump log's Gregory reference was not accepted")
    cfg = advice.get("config") or {}
    imbalance = cfg.get("imbalance", "energy_bot")
    heat_ratio = float(cfg.get("heat_ratio") or 1.0)
    applied = float(jump_log.get("somtp_dT_applied_mean", jump_log.get("somtp_dT")))
    dTS_exp = applied / heat_ratio

    t, N = _annual(columns, imbalance, "native")
    _, T = _annual(columns, "TS", "native")
    years = model_years(t, start_year)
    post = years >= jump_year
    settled = years >= jump_year + config.settle_years
    n_post, n_set = int(post.sum()), int(settled.sum())
    reasons: List[str] = []
    metrics: Dict[str, float] = {"dTS_expected": dTS_exp, "somtp_dT_applied": applied}

    def done(verdict):
        return CheckResult(verdict, reasons, jump_year, n_post, n_set, metrics)

    if n_post == 0:
        reasons.append(f"no complete model year since the jump at year {jump_year} "
                       f"in the trend data (latest is {int(years[-1])})")
        return done(Verdict.WAIT)
    if not (np.all(np.isfinite(N[post])) and np.all(np.isfinite(T[post]))):
        reasons.append("non-finite TS or N in the post-jump years: the run or its "
                       "trend output is broken")
        return done(Verdict.FAIL)

    pre = (years < jump_year) & (years >= jump_year - 10)
    sig_T = _noise(T[pre]) if pre.sum() >= 3 else _noise(T[post])
    sig_N = _noise(N[pre]) if pre.sum() >= 3 else _noise(N[post])
    # the no-jump first year: the pre-jump trend over the last baseline years
    last_pre = (years < jump_year) & (years >= jump_year - config.baseline_years)
    drift = 0.0
    if last_pre.sum() >= 3:
        drift = float(np.polyfit(t[last_pre], T[last_pre], 1)[0])  # K/yr, one year's worth
    # the jump's own effect relaxes during that first year (one box, tau):
    # its annual mean is tau*(1 - exp(-1/tau)) of the step
    tau = advice.get("tau_years")
    m = float(tau * (1.0 - np.exp(-1.0 / tau))) if tau and tau > 0 else 1.0

    # landed? (did the edited somtp reach the model?)
    moved = (float(T[post][0]) - TS_before - drift) / m
    metrics.update(TS_first=float(T[post][0]), TS_before=TS_before, drift_K_per_yr=drift,
                   dTS_first=moved, land_fraction=moved / dTS_exp if dTS_exp else float("nan"))
    if moved * dTS_exp > 0:
        metrics["heat_ratio_implied"] = applied / moved
    if (moved / dTS_exp < config.min_land_fraction
            and abs(dTS_exp - moved) > config.land_sigmas * sig_T):
        reasons.append(f"jump did not land: first post-jump TS moved {moved:+.2f} K "
                       f"(net of drift) vs {dTS_exp:+.2f} expected; was the edited "
                       f"docn.r the one the run read?")
        return done(Verdict.FAIL)

    if n_set == 0:
        reasons.append(f"jump landed; in the {config.settle_years}-yr adjustment "
                       f"period ({n_post} yr so far)")
        return done(Verdict.WAIT)

    last = (years < jump_year) & (years >= jump_year - config.baseline_years)
    bias = float(np.mean(N[last] - (T[last] - c0) / c1)) if last.any() else 0.0
    N_line = (T[settled] - c0) / c1
    dN = float(np.mean(N[settled] - N_line)) - bias
    tol_N = max(config.tol_N, 2.0 * sig_N / np.sqrt(n_set))
    metrics.update(N_obs=float(N[settled].mean()), N_line=float(N_line.mean()), dN=dN,
                   tol_N=tol_N, N_bias_before=bias, TS_obs=float(T[settled].mean()),
                   TS_expected=float(advice.get("TS_after") or np.nan))
    bad_N = abs(dN) > tol_N
    if bad_N:
        reasons.append(f"N is off the Gregory line by {dN:+.2f} W/m2 (tolerance "
                       f"{tol_N:.2f}): the run is not following the relation the "
                       f"jump was sized on")
    if n_set < config.min_years:
        if abs(dN) > config.early_fail_factor * tol_N:
            return done(Verdict.FAIL)
        reasons.append(f"{n_set} of {config.min_years} settled years so far")
        return done(Verdict.WAIT)
    return done(Verdict.FAIL if bad_N else Verdict.PASS)


@dataclass(frozen=True)
class ProbeCheckConfig:
    settle_years: int = 2
    #: settled years before a verdict other than WAIT
    min_years: int = 5
    #: after this many settled years an unreadable response is reported as
    #: PASS (lambda below what the probe can resolve) rather than WAIT forever
    max_wait_years: int = 12
    min_land_fraction: float = 0.3
    land_sigmas: float = 2.0
    #: significance of the imbalance change, in standard errors: n_se for a
    #: readable lambda, n_se_runaway for the runaway warning (stricter: the pre-probe
    #: mean rests on only recent_years years)
    n_se: float = 2.0
    n_se_runaway: float = 3.0
    #: the imbalance that decides the energy readout: energy_top, the planet's
    #: TOA balance (the quantity that defines convergence); energy_bot is
    #: read alongside for information only. Falls back to the advice's
    #: imbalance when the trend output lacks it.
    primary_imbalance: str = "energy_top"
    # ---- TS-trajectory readout (ts_relaxation_readout)
    #: pre-probe years the TS rate is fitted over (fewer than ts_pre_min_years
    #: usable -> no TS readout, the energy_bot readout stands alone)
    ts_pre_years: int = 20
    ts_pre_min_years: int = 8
    #: most recent settled post-probe years the post rate is fitted over (a
    #: straight line is only a good local fit of the relaxation while the
    #: window is short against tau)
    ts_post_max_years: int = 15
    #: significance for "TS rate changed" (readable) and for "TS rate rose after
    #: a warm probe" (runaway warning; stricter: a false alarm is costly)
    ts_n_se: float = 2.0
    ts_n_se_runaway: float = 3.0
    #: cap on the lag-1 autocorrelation used to inflate standard errors
    ts_rho_max: float = 0.5
    #: heat capacity (W yr m-2 K-1) behind lambda = C/tau; None = take it from
    #: the advice (C_eff_W_yr_m2_K, else heat.C_ocean)
    heat_capacity: Optional[float] = None

    def __post_init__(self) -> None:
        for name in ("settle_years", "min_years", "max_wait_years", "ts_pre_years",
                     "ts_pre_min_years", "ts_post_max_years"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")


def _line(x: np.ndarray, y: np.ndarray) -> dict:
    """OLS line with the centred-x quantities the readout needs."""
    xm, ym = float(x.mean()), float(y.mean())
    sxx = float(((x - xm) ** 2).sum())
    slope = float(((x - xm) * (y - ym)).sum() / sxx)
    res = y - ym - slope * (x - xm)
    return {"n": int(x.size), "xm": xm, "ym": ym, "sxx": sxx, "slope": slope,
            "res": res, "ssr": float((res ** 2).sum())}


def _tau_eq(Tp: float, Tq: float, rp: float, rq: float):
    """One-box relaxation from two (T, dT/dt) points: tau and the equilibrium."""
    D = rp - rq
    tau = (Tq - Tp) / D
    return tau, Tp + rp * tau


def ts_relaxation_readout(years, T, jump_year: int, settle_years: int = 2,
                          config: ProbeCheckConfig = ProbeCheckConfig(),
                          probe_sign: float = 1.0, since_year: Optional[int] = None,
                          heat_capacity: Optional[float] = None) -> Optional[dict]:
    """Read a probe from the surface-temperature trajectory (annual TS).

    energy_bot is a noisy gauge on a hot slab aquaplanet (annual sigma ~8 W/m2
    on D2, so five settled years resolve only lambda >~ 2 W/m2/K) while annual
    TS scatters by ~1 K, and its drift is measurable both before and after the
    probe. The one-box relaxation ``dTS/dt = (TS_eq - TS)/tau`` gives two
    (rate, temperature) points on the same line:

        r_pre  at T_pre   (OLS over the last ``ts_pre_years`` pre-probe years)
        r_post at T_post  (OLS over the settled post-probe years, at most
                           ``ts_post_max_years``)
        tau    = (T_post - T_pre) / (r_pre - r_post)
        TS_eq  = T_pre + r_pre * tau

    Why two local slopes rather than a joint exponential fit of the whole
    series: the post-probe series is a transient, but over a window short
    against tau the OLS slope is the derivative at the window's mid-time and
    the window mean is the temperature there (error ~ (L/tau)^2/24 of the
    rate, ~1.5 % for L=10 yr, tau=25 yr), so each (T, r) pair is unbiased.
    The estimators are then (T, r) pairs with Gaussian errors, which
    propagate through the two equations by the delta method; a joint
    nonlinear fit of (TS_eq, tau, landed step) from a 20-yr pre window that
    is nearly straight is badly conditioned and its errors are not Gaussian.
    The information used is the same: the *change* in rate across a known
    change in temperature.

    Uncertainties. Slope standard errors use the pooled residual scatter of
    both fits (the post residuals carry a little curvature, which only
    inflates them, conservatively) and are inflated by (1+rho)/(1-rho)
    for the lag-1 autocorrelation rho of the *pre-probe* residuals (more
    degrees of freedom; clipped to [0, ``ts_rho_max``]: a 20-yr estimate of rho
    is itself noisy, and the cap keeps one lucky run of years from inflating
    the error without bound). Means: sigma_eff/sqrt(n). The four inputs are
    independent (centred OLS gives uncorrelated mean and slope). tau and
    TS_eq errors are first order and trustworthy only when the rate change is
    well resolved; below ``ts_n_se`` sigma they are not reported, only a
    lower bound on tau. The uncertainty of the heat capacity is not included.

    Returns None when the pre-probe window is too short, else a dict with the
    ``ts_`` metrics (status 'runaway', 'readable', 'falling_back' or
    'unresolved'; 'side' 'bracketed' / 'above' where known).
    """
    years = np.asarray(years)
    T = np.asarray(T, dtype=float)
    pre = (years < jump_year) & (years >= jump_year - config.ts_pre_years)
    if since_year is not None:
        pre &= years >= since_year
    post = years >= jump_year + settle_years
    post_idx = np.flatnonzero(post)[-config.ts_post_max_years:]
    post = np.zeros(years.size, dtype=bool)
    post[post_idx] = True
    if pre.sum() < config.ts_pre_min_years or post.sum() < 3:
        return None
    if not (np.all(np.isfinite(T[pre])) and np.all(np.isfinite(T[post]))):
        return None
    a = _line(years[pre].astype(float), T[pre])
    b = _line(years[post].astype(float), T[post])
    dof = (a["n"] - 2) + (b["n"] - 2)
    sigma = float(np.sqrt((a["ssr"] + b["ssr"]) / dof))
    e = a["res"]
    rho = float((e[:-1] * e[1:]).sum() / (e ** 2).sum()) if e.size > 2 else 0.0
    rho = min(max(rho, 0.0), config.ts_rho_max)
    s_eff = sigma * np.sqrt((1.0 + rho) / (1.0 - rho))
    se = {"Tp": s_eff / np.sqrt(a["n"]), "Tq": s_eff / np.sqrt(b["n"]),
          "rp": s_eff / np.sqrt(a["sxx"]), "rq": s_eff / np.sqrt(b["sxx"])}
    Tp, Tq, rp, rq = a["ym"], b["ym"], a["slope"], b["slope"]
    s = 1.0 if probe_sign >= 0 else -1.0
    D = rp - rq
    se_D = float(np.hypot(se["rp"], se["rq"]))
    z = s * D / se_D                                  # >0: the probe slowed the drift
    out = {"ts_r_pre": rp, "ts_r_pre_se": se["rp"], "ts_r_post": rq,
           "ts_r_post_se": se["rq"], "ts_T_pre_fit": Tp, "ts_T_post_fit": Tq,
           "ts_sigma": sigma, "ts_rho": rho, "ts_n_pre": float(a["n"]),
           "ts_n_post": float(b["n"]), "ts_D": D, "ts_D_se": se_D, "ts_D_sigmas": z}

    if z < -config.ts_n_se_runaway:
        out["ts_status"] = "runaway"
        return out
    falling = s * rq < -config.ts_n_se * se["rq"]    # TS reverses against the probe
    if falling:
        out["ts_side"] = "bracketed"
    if z >= config.ts_n_se and s * (Tq - Tp) > 0:
        tau, Teq = _tau_eq(Tp, Tq, rp, rq)
        # first-order (delta-method) errors from numerical gradients
        base = np.array([Tp, Tq, rp, rq])
        keys = ("Tp", "Tq", "rp", "rq")
        var_tau = var_eq = 0.0
        for i, k in enumerate(keys):
            h = 1e-6 * max(1.0, abs(base[i]))
            up, dn = base.copy(), base.copy()
            up[i] += h
            dn[i] -= h
            tu, eu = _tau_eq(*up)
            td, ed = _tau_eq(*dn)
            var_tau += ((tu - td) / (2 * h) * se[k]) ** 2
            var_eq += ((eu - ed) / (2 * h) * se[k]) ** 2
        out.update(ts_status="readable", ts_tau=tau, ts_tau_se=float(np.sqrt(var_tau)),
                   ts_eq=Teq, ts_eq_se=float(np.sqrt(var_eq)),
                   ts_side="above" if s * (Teq - Tq) > 0 else "bracketed")
        if heat_capacity:
            out["ts_lambda"] = heat_capacity / tau
            out["ts_lambda_se"] = heat_capacity / tau ** 2 * out["ts_tau_se"]
        return out
    if falling:
        out["ts_status"] = "falling_back"
        return out
    out["ts_status"] = "unresolved"
    denom = s * D + config.ts_n_se * se_D          # rate change at its upper edge
    if denom > 0 and s * (Tq - Tp) > 0:
        out["ts_tau_lower"] = float(abs(Tq - Tp) / denom)
        if heat_capacity:
            out["ts_lambda_upper"] = heat_capacity / out["ts_tau_lower"]
    return out


def _heat_capacity(advice: dict, config: ProbeCheckConfig) -> Optional[float]:
    if config.heat_capacity:
        return float(config.heat_capacity)
    for v in (advice.get("C_eff_W_yr_m2_K"), (advice.get("heat") or {}).get("C_ocean")):
        if v is not None and np.isfinite(v) and v > 0:
            return float(v)
    return None


def check_ocean_probe(columns: Dict[str, np.ndarray], jump_log: dict,
                      start_year: int = 1,
                      config: ProbeCheckConfig = ProbeCheckConfig()) -> CheckResult:
    """Read a probe's response: lambda, the implied equilibrium, and the side.

    Two independent readouts, each reported, then combined:

    *Energy readout* (``primary_imbalance``, default energy_top — the TOA
    balance that defines convergence; energy_bot is read the same way and
    reported for information): compares the settled post-probe years with the
    ``recent_years`` before the probe (the probe's lever arm), the standard
    error from the interannual scatter of a long detrended pre-probe window:

    * lambda = -(N_post - N_pre) / (TS_post - TS_pre), with its standard error;
    * runaway — the imbalance *grew* in the direction of the probe (a warm probe
      that increases energy_bot): no restoring feedback at the new state;
    * readable — lambda significantly > 0: ``TS_eq = TS_post + N_post/lambda``
      and ``side`` (N_post > 0 after a warm probe: the equilibrium is hotter
      still; < 0: the probe overshot and the equilibrium is bracketed between
      the pre- and post-probe states).

    *TS readout* (``ts_relaxation_readout``): the surface-temperature drift
    before vs after the probe -> tau, TS_eq, side, and lambda = C/tau when the
    advice carries a heat capacity. Far quieter than energy_bot on noisy runs.

    Combination (conservative, explicit):

    * FAIL — the probe never landed or the data are non-finite;
    * warning "runaway greenhouse suspected" (verdict WAIT, no recommendation
      either way) — the energy readout grew in the probe's direction by
      ``n_se_runaway`` (3) sigma, or the TS drift *rose* after a warm probe by
      ``ts_n_se_runaway`` (3) sigma;
    * PASS — at least one readout is significant (the energy one at ``n_se``,
      the TS one at ``ts_n_se`` sigma, or TS reversing against the probe) and
      no FAIL rule fires and they do not contradict. Two readouts contradict
      when both are readable and put the equilibrium on different sides of
      the probe level (hotter still vs bracketed): the verdict is then WAIT;
    * WAIT — fewer than ``min_years`` settled years, neither readout
      significant, or a contradiction. After ``max_wait_years`` settled years
      with neither readable: PASS, with upper bounds on |lambda| / a lower
      bound on tau (the response is below what the probe resolves).
    """
    advice = jump_log.get("advice") or {}
    cfg = advice.get("config") or {}
    imbalance = cfg.get("imbalance", "energy_bot")
    k = int(cfg.get("recent_years") or 5)
    jump_year = int(jump_log["jump_model_year"])
    applied = float(jump_log.get("somtp_dT_applied_mean", jump_log.get("somtp_dT")))
    hr = float(cfg.get("heat_ratio") or 1.0)
    expected = applied / hr

    t, N = _annual(columns, imbalance, "native")
    _, T = _annual(columns, "TS", "native")
    years = model_years(t, start_year)
    post = years >= jump_year
    settled = years >= jump_year + config.settle_years
    pre = (years < jump_year) & (years >= jump_year - k)
    n_post, n_set = int(post.sum()), int(settled.sum())
    reasons: List[str] = []
    metrics: Dict[str, float] = {"somtp_dT_applied": applied, "dTS_expected": expected}

    warnings: List[str] = []

    def done(verdict):
        return CheckResult(verdict, reasons, jump_year, n_post, n_set, metrics, warnings)

    if n_post == 0:
        reasons.append(f"no complete model year since the probe at year {jump_year}")
        return done(Verdict.WAIT)
    if not (np.all(np.isfinite(N[post | pre])) and np.all(np.isfinite(T[post | pre]))):
        reasons.append("non-finite TS or N around the probe: the run or its trend "
                       "output is broken")
        return done(Verdict.FAIL)
    if pre.sum() < 3:
        raise ValueError(f"only {int(pre.sum())} pre-probe years in the trend data")

    sig_T = _noise(T[pre | (years >= jump_year - 10) & (years < jump_year)])
    drift = float(np.polyfit(t[pre], T[pre], 1)[0])
    TS_before = float(advice.get("TS_now", T[pre][-1]))
    moved = float(T[post][0]) - TS_before - drift
    metrics.update(TS_first=float(T[post][0]), dTS_first=moved,
                   land_fraction=moved / expected if expected else float("nan"))
    if moved * applied > 0:
        metrics["heat_ratio_implied"] = applied / moved
    if (moved / expected < config.min_land_fraction
            and abs(expected - moved) > config.land_sigmas * sig_T):
        reasons.append(f"probe did not land: first post-probe TS moved {moved:+.2f} K "
                       f"(net of drift) vs {expected:+.2f}; was the edited docn.r the "
                       f"one the run read?")
        return done(Verdict.FAIL)
    if n_set < config.min_years:
        reasons.append(f"probe landed; {n_set} of {config.min_years} settled years "
                       f"so far")
        return done(Verdict.WAIT)

    # ---- energy readouts: energy_top decides, energy_bot is information
    # interannual noise of the imbalance comes from detrended residuals over a
    # long pre-probe window (and the settled post years), not from the scatter
    # of the recent_years averaged: a 5-sample variance is itself so noisy that
    # a 2-3 sigma rule on it false-alarms several times too often
    noise_pre = (years < jump_year) & (years >= jump_year - config.ts_pre_years)
    if advice.get("since_year") is not None:
        noise_pre &= years >= int(advice["since_year"])
    primary = (config.primary_imbalance
               if _have(columns, config.primary_imbalance, "native") else imbalance)
    readouts = {}
    for var in dict.fromkeys((primary, "energy_bot", "energy_top")):
        if not _have(columns, var, "native"):
            continue
        _, Nv = _annual(columns, var, "native")
        if not np.all(np.isfinite(Nv[post | pre])):
            reasons.append(f"non-finite {var} around the probe: the run or its trend "
                           f"output is broken")
            return done(Verdict.FAIL)
        readouts[var] = _energy_readout(Nv, T, pre, settled, config.n_se,
                                        config.n_se_runaway, noise_pre)
    er = readouts[primary]
    lam, lam_se, dN, se = er["lambda"], er["lambda_se"], er["dN"], er["se_dN"]
    metrics.update(N_pre=er["N_pre"], N_post=er["N_post"], dN=dN, se_dN=se,
                   TS_pre=er["TS_pre"], TS_post=er["TS_post"], dTS=er["dTS"],
                   **{"lambda": lam, "lambda_se": lam_se})
    metrics["primary_is_top"] = 1.0 if primary == "energy_top" else 0.0
    s = 1.0 if er["dTS"] > 0 else -1.0
    n_state, n_side = er["state"], er["side"]
    if n_state == "readable":
        metrics["TS_eq"] = er["TS_eq"]
        metrics["bracketed"] = 1.0 if n_side == "bracketed" else 0.0
    reasons.append(_energy_sentence(primary, er))
    for var, r in readouts.items():
        if var == primary:
            continue
        tag = "bot" if var == "energy_bot" else "top"
        metrics.update({f"{tag}_lambda": r["lambda"], f"{tag}_lambda_se": r["lambda_se"],
                        f"{tag}_dN": r["dN"], f"{tag}_N_post": r["N_post"]})
        reasons.append(f"(information) {_energy_sentence(var, r)}")

    # ---- TS readout
    C = _heat_capacity(advice, config)
    if C:
        metrics["heat_capacity"] = C
    ts = ts_relaxation_readout(years, T, jump_year, config.settle_years, config,
                               probe_sign=s, since_year=advice.get("since_year"),
                               heat_capacity=C)
    ts_state = ts["ts_status"] if ts else "none"
    if ts:
        metrics.update({k_: v for k_, v in ts.items() if isinstance(v, float)})
        if "ts_side" in ts:
            metrics["ts_bracketed"] = 1.0 if ts["ts_side"] == "bracketed" else 0.0
        reasons.append(_ts_sentence(ts, C))

    # A probe that is pushed further away instead of pulled back (the
    # imbalance grows in the probe's direction, or TS warms faster than before
    # a warm probe) has no restoring feedback at the new state. That is
    # reported as a warning, not a verdict: whether to keep running or roll
    # back is the user's decision.
    if n_state == "runaway":
        warnings.append(f"runaway greenhouse suspected: {primary} rose by "
                        f"{dN:+.2f} ± {se:.2f} W/m2 after a {er['dTS']:+.2f} K probe "
                        f"(>= {config.n_se_runaway:g} sigma) — no restoring feedback seen")
    if ts_state == "runaway":
        warnings.append(f"runaway greenhouse suspected: TS drift rose after the probe "
                        f"({ts['ts_r_pre']:+.3f} -> {ts['ts_r_post']:+.3f} K/yr, "
                        f"{-ts['ts_D_sigmas']:.1f} sigma) — TS accelerating away")
    if warnings:
        metrics["runaway_suspected"] = 1.0
        reasons.append("runaway greenhouse suspected (see warnings): no PASS while "
                       "it stands; keep watching or roll back — your call")
        return done(Verdict.WAIT)
    n_ok = n_state == "readable"
    ts_ok = ts_state in ("readable", "falling_back")
    ts_side = ts.get("ts_side") if ts else None
    if n_ok and ts_ok and ts_side and n_side != ts_side:
        reasons.append(f"the readouts contradict: {primary} puts the equilibrium "
                       f"'{n_side}', the TS drift '{ts_side}' — keep running")
        return done(Verdict.WAIT)
    if n_ok or ts_ok:
        return done(Verdict.PASS)
    if n_set >= config.max_wait_years:
        msg = (f"response within the noise after {n_set} settled years: "
               f"|lambda| < ~{config.n_se * lam_se:.2f} W/m2/K at this state (energy)")
        if ts and "ts_tau_lower" in ts:
            msg += f"; tau > {ts['ts_tau_lower']:.0f} yr (TS)"
        reasons.append(msg)
        return done(Verdict.PASS)
    reasons.append(f"{primary} changed {dN:+.2f} ± {se:.2f} W/m2 so far: not yet "
                   f"readable (lambda {lam:+.2f} ± {lam_se:.2f})")
    return done(Verdict.WAIT)


def _detrended_ssr(t, y):
    if y.size < 3:
        return 0.0, 0
    r = y - np.polyval(np.polyfit(t, y, 1), t)
    return float((r ** 2).sum()), y.size - 2


def _energy_readout(N, T, pre, settled, n_se: float, n_se_runaway: float,
                    noise_pre=None) -> dict:
    """Lever-arm readout of one imbalance series: settled post-probe years
    against the pre-probe years. state 'runaway' (the imbalance grew in the
    probe's direction), 'readable' (lambda significantly > 0) or 'noise'.
    The standard error uses the interannual scatter pooled from detrended
    residuals over ``noise_pre`` (a longer pre-probe window) and the settled
    years, when that window is longer than ``pre``."""
    Na, Nb = N[pre], N[settled]
    TS_pre, TS_post = float(T[pre].mean()), float(T[settled].mean())
    dTS = TS_post - TS_pre
    dN = float(Nb.mean() - Na.mean())
    idx = np.arange(N.size, dtype=float)
    if noise_pre is not None and noise_pre.sum() > pre.sum():
        a, da = _detrended_ssr(idx[noise_pre], N[noise_pre])
        b, db = _detrended_ssr(idx[settled], Nb)
        sig = float(np.sqrt((a + b) / (da + db))) if da + db > 0 else float("nan")
        se = sig * float(np.sqrt(1.0 / Na.size + 1.0 / Nb.size))
    else:
        se = float(np.sqrt(Na.var(ddof=1) / Na.size + Nb.var(ddof=1) / Nb.size))
    s = 1.0 if dTS > 0 else -1.0
    out = {"N_pre": float(Na.mean()), "N_post": float(Nb.mean()), "dN": dN, "se_dN": se,
           "TS_pre": TS_pre, "TS_post": TS_post, "dTS": dTS,
           "lambda": -dN / dTS if dTS else float("nan"),
           "lambda_se": se / abs(dTS) if dTS else float("nan"),
           "state": "noise", "side": None}
    if s * dN > n_se_runaway * se:
        out["state"] = "runaway"
    elif -s * dN > n_se * se:
        out["state"] = "readable"
        out["TS_eq"] = TS_post + out["N_post"] / out["lambda"]
        out["side"] = "above" if out["N_post"] * s > 0 else "bracketed"
    return out


def _energy_sentence(var: str, r: dict) -> str:
    lam = f"lambda {r['lambda']:.2f} ± {r['lambda_se']:.2f} W/m2/K"
    if r["state"] == "runaway":
        return (f"{var} rose by {r['dN']:+.2f} ± {r['se_dN']:.2f} W/m2 after a "
                f"{r['dTS']:+.2f} K probe: no restoring feedback at this state")
    if r["state"] == "readable":
        if r["side"] == "above":
            side = "hotter still" if r["dTS"] > 0 else "cooler still"
            return (f"{var}: {lam}, still {r['N_post']:+.2f} W/m2: the equilibrium is "
                    f"{side}, ~{r['TS_eq']:.1f} K — another probe or a Gregory jump "
                    f"can follow")
        return (f"{var}: {lam}, now {r['N_post']:+.2f} W/m2: the probe overshot — the "
                f"equilibrium (~{r['TS_eq']:.1f} K) lies between {r['TS_pre']:.1f} and "
                f"{r['TS_post']:.1f} K")
    return f"{var} changed {r['dN']:+.2f} ± {r['se_dN']:.2f} W/m2 ({lam}): within the noise"

def _ts_sentence(ts: dict, C: Optional[float]) -> str:
    base = (f"TS drift {ts['ts_r_pre']:+.3f} ± {ts['ts_r_pre_se']:.3f} K/yr before, "
            f"{ts['ts_r_post']:+.3f} ± {ts['ts_r_post_se']:.3f} after")
    st = ts["ts_status"]
    if st == "readable":
        lam = (f", lambda {ts['ts_lambda']:.2f} ± {ts['ts_lambda_se']:.2f} W/m2/K"
               if "ts_lambda" in ts else "")
        where = ("equilibrium above the probe level" if ts["ts_side"] == "above"
                 else "TS falls back: equilibrium below the probe level (bracketed)")
        return (f"{base}: tau {ts['ts_tau']:.0f} ± {ts['ts_tau_se']:.0f} yr, TS_eq "
                f"{ts['ts_eq']:.1f} ± {ts['ts_eq_se']:.1f} K{lam}; {where}")
    if st == "falling_back":
        return f"{base}: TS is falling back — equilibrium below the probe level (bracketed)"
    if st == "runaway":
        return base
    extra = (f"; tau > {ts['ts_tau_lower']:.0f} yr" if "ts_tau_lower" in ts else "")
    return f"{base}: change not yet resolved ({ts['ts_D_sigmas']:+.1f} sigma){extra}"


def check_any(columns, jump_log: dict, start_year: int = 1,
              settle_years: Optional[int] = None,
              min_years: Optional[int] = None) -> CheckResult:
    """Dispatch on the jump log's plugin (aqua_ice when unrecorded)."""
    kw = {k: v for k, v in (("settle_years", settle_years), ("min_years", min_years))
          if v is not None}
    if jump_log.get("plugin") == "som_ocean":
        return check_ocean_jump(columns, jump_log, start_year, OceanCheckConfig(**kw))
    return check_jump(columns, jump_log, start_year, CheckConfig(**kw))
