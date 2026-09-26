"""Post-jump consistency check: did the run continue along the predicted path?

Pure computation. Implements the design of ``consistency.py`` for the aqua_ice
jump. Inputs are the case's exocam-trend columns after the jump and the jump
log (the ``<cice.r>.accel.json`` sidecar ``restart.apply_ice_jump`` writes),
which holds the advice the jump was sized from.

The reference is the advice's own fitted relations, evaluated at the ice the
run actually has each year — not a single number:

* landed:   first post-jump annual hi ~ hi_before x factor (else the edit
            never reached the model: wrong file, stale rpointer, rollback);
* on-law:   N_obs ~ a + b/hi_obs   (the conduction law the jump relied on);
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
    factor = float(jump_log["ice_factor"])
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
    N_law = a + b / h[settled]
    dN = float(np.mean(N[settled] - N_law)) - bias_N
    tol_N = max(config.tol_N, 2.0 * sig_N / np.sqrt(n_set))
    metrics.update(N_obs=float(N[settled].mean()), N_law=float(N_law.mean()),
                   dN=dN, tol_N=tol_N, N_bias_before=bias_N)

    # TS on the Gregory relation?
    dT = tol_T = None
    ts = (advice.get("temperatures") or {}).get("TS") or {}
    lin = (ts.get("fits") or {}).get("linear")
    if lin and f"TS_native" in columns:
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
    if dT is None:
        reasons.append("no TS reference in the jump log: TS not checked")

    melting = "hi_trend_m_per_yr" in metrics and metrics["hi_trend_m_per_yr"] < 0
    if n_set < config.min_years:
        k = config.early_fail_factor
        if abs(dN) > k * tol_N or (dT is not None and abs(dT) > k * tol_T):
            return done(Verdict.FAIL)
        reasons.append(f"{n_set} of {config.min_years} settled years so far")
        return done(Verdict.WAIT)
    if bad_N or bad_T or melting:
        return done(Verdict.FAIL)
    return done(Verdict.PASS)
