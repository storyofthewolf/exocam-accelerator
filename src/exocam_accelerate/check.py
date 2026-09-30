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

from .advise import ICE_VOLUME_VAR, _annual, model_years


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
    #: significance of the energy_bot change, in standard errors
    n_se: float = 2.0

    def __post_init__(self) -> None:
        for name in ("settle_years", "min_years", "max_wait_years"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")


def check_ocean_probe(columns: Dict[str, np.ndarray], jump_log: dict,
                      start_year: int = 1,
                      config: ProbeCheckConfig = ProbeCheckConfig()) -> CheckResult:
    """Read a probe's response: lambda, the implied equilibrium, and the side.

    Compares the settled post-probe years with the ``recent_years`` before
    the probe (the probe's lever arm):

    * lambda = -(N_post - N_pre) / (TS_post - TS_pre), with its standard error;
    * FAIL — the imbalance *grew* in the direction of the probe (a warm probe
      that increases energy_bot): no restoring feedback at the new state, the
      runaway signature; also a probe that never landed, or non-finite data;
    * PASS — lambda significantly > 0: ``TS_eq = TS_post + N_post/lambda`` and
      ``side`` (N_post > 0 after a warm probe: the equilibrium is hotter still;
      < 0: the probe overshot and the equilibrium is bracketed between the
      pre- and post-probe states);
    * WAIT — fewer than ``min_years`` settled years, or a response still inside
      the noise (until ``max_wait_years``, then PASS with lambda below what the
      probe resolves).
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

    def done(verdict):
        return CheckResult(verdict, reasons, jump_year, n_post, n_set, metrics)

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

    Na, Nb = N[pre], N[settled]
    dTS = float(T[settled].mean() - T[pre].mean())
    dN = float(Nb.mean() - Na.mean())
    se = float(np.sqrt(Na.var(ddof=1) / Na.size + Nb.var(ddof=1) / Nb.size))
    s = 1.0 if dTS > 0 else -1.0
    lam = -dN / dTS if dTS else float("nan")
    lam_se = se / abs(dTS) if dTS else float("nan")
    N_post, TS_post = float(Nb.mean()), float(T[settled].mean())
    metrics.update(N_pre=float(Na.mean()), N_post=N_post, dN=dN, se_dN=se,
                   TS_pre=float(T[pre].mean()), TS_post=TS_post, dTS=dTS,
                   **{"lambda": lam, "lambda_se": lam_se})

    if s * dN > config.n_se * se:
        reasons.append(f"{imbalance} rose by {dN:+.2f} ± {se:.2f} W/m2 after a "
                       f"{dTS:+.2f} K probe: no restoring feedback at this state "
                       f"(runaway-like) — roll back")
        return done(Verdict.FAIL)
    if -s * dN > config.n_se * se:
        TS_eq = TS_post + N_post / lam
        metrics["TS_eq"] = TS_eq
        if N_post * s > 0:
            side = "hotter still" if s > 0 else "cooler still"
            metrics["bracketed"] = 0.0
            reasons.append(f"lambda {lam:.2f} ± {lam_se:.2f} W/m2/K; {imbalance} still "
                           f"{N_post:+.2f}: the equilibrium is {side}, ~{TS_eq:.1f} K "
                           f"— another probe or a Gregory jump can follow")
        else:
            metrics["bracketed"] = 1.0
            reasons.append(f"lambda {lam:.2f} ± {lam_se:.2f} W/m2/K; {imbalance} "
                           f"{N_post:+.2f}: the probe overshot — the equilibrium "
                           f"(~{TS_eq:.1f} K) lies between {metrics['TS_pre']:.1f} and "
                           f"{TS_post:.1f} K")
        return done(Verdict.PASS)
    if n_set >= config.max_wait_years:
        reasons.append(f"response within the noise after {n_set} settled years: "
                       f"|lambda| < ~{config.n_se * lam_se:.2f} W/m2/K at this state")
        return done(Verdict.PASS)
    reasons.append(f"{imbalance} changed {dN:+.2f} ± {se:.2f} W/m2 so far: not yet "
                   f"readable (lambda {lam:+.2f} ± {lam_se:.2f})")
    return done(Verdict.WAIT)


def check_any(columns, jump_log: dict, start_year: int = 1,
              settle_years: Optional[int] = None,
              min_years: Optional[int] = None) -> CheckResult:
    """Dispatch on the jump log's plugin (aqua_ice when unrecorded)."""
    kw = {k: v for k, v in (("settle_years", settle_years), ("min_years", min_years))
          if v is not None}
    if jump_log.get("plugin") == "som_ocean":
        return check_ocean_jump(columns, jump_log, start_year, OceanCheckConfig(**kw))
    return check_jump(columns, jump_log, start_year, CheckConfig(**kw))
