"""Phase-space extrapolation: a state variable against the energy imbalance N.

Pure computation (numpy only). Motivation: docs/phase-space-extrapolation.md.
Fitting a state X against the net TOA imbalance N and evaluating at a chosen
N is far better conditioned than extrapolating X against time, because the
target abscissa is chosen (or pinned by physics) rather than inferred from a
timescale.

Three fit forms:

  linear      X = c0 + c1*N                     Gregory-style; the reference
                                                for temperatures (TS, Tsfc).
  saturating  X = X0 - A*(1 - exp(N/N0))        curved temperature response
                                                (bake-off "phase-nl").
  hyperbolic  N = a + b/X   <=>  X = b/(N - a)  conduction-limited (Stefan)
                                                ice growth: the deficit N is
                                                conducted through ice of
                                                thickness X, so N*X ~ const.
                                                X diverges as N -> a (~0):
                                                thick ice never "equilibrates";
                                                a jump targets an N, not N=0.

The hyperbolic law was found in the 15 grp3 cold cases (N*hi constant to a
few percent over yr 30-150, r^2 >= 0.98 in 12/15; docs "Ice growth is
Stefan-limited").

Gates refuse, never shrink:
  * too few points / no variation in the window,
  * X not slaved to N (|corr| small in the fitted coordinates),
  * target N too far beyond the sampled N-range,
  * linear vs saturating endpoints disagree (temperatures: a kink),
  * hyperbolic target at or past the asymptote a (X would be infinite).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

# N0 search range (W/m2) for the saturating form. The upper end is effectively
# linear over any realistic imbalance range; the lower end forbids curvature
# on scales finer than interannual noise in N.
_N0_MIN = 0.5
_N0_MAX = 1.0e4

FORMS = ("linear", "saturating", "hyperbolic")


@dataclass(frozen=True)
class PhaseFit:
    """One fitted relation X(N)."""

    form: str                    # one of FORMS
    params: Tuple[float, ...]    # linear (c0, c1); saturating (X0, A, N0);
                                 # hyperbolic (a, b)
    rmse: float                  # in X for linear/saturating, in N for hyperbolic
    corr: float                  # corr in the fitted coordinates
    at_bound: bool = False       # saturating N0 pinned at its search bound

    def predict(self, N) -> np.ndarray:
        N = np.asarray(N, dtype=float)
        if self.form == "linear":
            c0, c1 = self.params
            return c0 + c1 * N
        if self.form == "saturating":
            X0, A, N0 = self.params
            return X0 - A * (1.0 - np.exp(N / N0))
        a, b = self.params
        with np.errstate(divide="ignore"):
            X = b / (N - a)
        return np.where(X > 0, X, np.inf)


def fit_linear(N, X) -> PhaseFit:
    N = np.asarray(N, dtype=float)
    X = np.asarray(X, dtype=float)
    c1, c0 = np.polyfit(N, X, 1)
    resid = X - (c0 + c1 * N)
    return PhaseFit("linear", (float(c0), float(c1)),
                    float(np.sqrt(np.mean(resid**2))),
                    float(np.corrcoef(N, X)[0, 1]))


def _saturating_lsq(N, X, N0):
    """Exact (X0, A) and SSE for fixed N0."""
    g = 1.0 - np.exp(N / N0)
    G = np.column_stack([np.ones_like(N), -g])
    coef, *_ = np.linalg.lstsq(G, X, rcond=None)
    resid = X - G @ coef
    return float(coef[0]), float(coef[1]), float(resid @ resid)


def fit_saturating(N, X, n_grid: int = 241) -> PhaseFit:
    """Fit X = X0 - A*(1 - exp(N/N0)): log-grid + golden refinement over N0.

    Linear in (X0, A) for fixed N0, so each trial is an exact least-squares
    solve — no scipy, no initial-guess sensitivity.
    """
    N = np.asarray(N, dtype=float)
    X = np.asarray(X, dtype=float)
    grid = np.logspace(np.log10(_N0_MIN), np.log10(_N0_MAX), n_grid)
    sse = np.array([_saturating_lsq(N, X, n0)[2] for n0 in grid])
    i = int(np.argmin(sse))

    lo = np.log(grid[max(i - 1, 0)])
    hi = np.log(grid[min(i + 1, n_grid - 1)])
    phi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc = _saturating_lsq(N, X, np.exp(c))[2]
    fd = _saturating_lsq(N, X, np.exp(d))[2]
    for _ in range(40):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = _saturating_lsq(N, X, np.exp(c))[2]
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = _saturating_lsq(N, X, np.exp(d))[2]
    n0 = float(np.exp((a + b) / 2.0))
    X0, A, s = _saturating_lsq(N, X, n0)
    return PhaseFit("saturating", (X0, A, n0), float(np.sqrt(s / N.size)),
                    float(np.corrcoef(N, X)[0, 1]),
                    at_bound=i == 0 or i == n_grid - 1)


def fit_hyperbolic(N, X) -> PhaseFit:
    """Fit N = a + b/X (linear regression of N on 1/X). Requires X > 0."""
    N = np.asarray(N, dtype=float)
    X = np.asarray(X, dtype=float)
    if np.any(X <= 0):
        raise ValueError("hyperbolic fit needs X > 0 throughout the window")
    inv = 1.0 / X
    b, a = np.polyfit(inv, N, 1)
    resid = N - (a + b * inv)
    return PhaseFit("hyperbolic", (float(a), float(b)),
                    float(np.sqrt(np.mean(resid**2))),
                    float(np.corrcoef(inv, N)[0, 1]))


@dataclass(frozen=True)
class PhaseGateConfig:
    """Thresholds for trusting a phase-space extrapolation.

    Defaults from the offline hindcast on the 15 grp3 cold cases
    (docs/phase-space-extrapolation.md, "Advisor hindcast").
    """

    min_points: int = 8
    #: |corr| in the fitted coordinates (N vs X, or N vs 1/X)
    min_abs_corr: float = 0.9
    #: (N_target - nearest sampled N) / (sampled N range), when beyond it
    max_extrapolation_ratio: float = 5.0
    #: |pred_linear - pred_saturating| / |step| (temperature forms only)
    max_relative_disagreement: float = 1.0

    def __post_init__(self) -> None:
        if self.max_extrapolation_ratio <= 0:
            raise ValueError(f"max_extrapolation_ratio must be positive, got "
                             f"{self.max_extrapolation_ratio!r}")
        if self.min_points < 2:
            raise ValueError(f"min_points must be at least 2, got {self.min_points!r}")
        if not (0.0 <= self.min_abs_corr <= 1.0):
            raise ValueError(f"min_abs_corr must satisfy 0 <= min_abs_corr <= 1, got "
                             f"{self.min_abs_corr!r}")


@dataclass(frozen=True)
class PhaseExtrapolation:
    """X extrapolated in phase space to a target imbalance, or a refusal."""

    variable: str
    form: str
    accepted: bool
    reasons: List[str]
    N_target: float
    X_now: float
    prediction: float = float("nan")      # NaN when refused
    fits: dict = field(default_factory=dict)          # form -> PhaseFit
    predictions: dict = field(default_factory=dict)   # form -> X(N_target)
    extrapolation_ratio: float = float("nan")

    @property
    def step(self) -> float:
        return self.prediction - self.X_now

    def __bool__(self) -> bool:
        return self.accepted


def _extrapolation_ratio(N, N_target) -> float:
    span = float(np.ptp(N))
    if N.min() <= N_target <= N.max():
        return 0.0
    return float(np.min(np.abs(N - N_target))) / span


def extrapolate(
    variable: str,
    N,
    X,
    N_target: float,
    X_now: float,
    form: str = "linear",
    config: PhaseGateConfig = PhaseGateConfig(),
) -> PhaseExtrapolation:
    """Fit X(N) over a window and evaluate at ``N_target``, behind the gate.

    ``N``/``X`` are (smoothed, annual-mean) window samples; ``X_now`` is the
    current state the step is measured from (normally the latest native annual
    mean, not the lagged running mean).

    ``form="linear"`` or ``"saturating"`` fits both temperature forms, returns
    the requested one, and refuses when they disagree. ``form="hyperbolic"``
    is the ice-reservoir form and refuses a target at or past its asymptote.
    """
    if form not in FORMS:
        raise ValueError(f"unknown form {form!r}; expected one of {FORMS}")
    N = np.asarray(N, dtype=float)
    X = np.asarray(X, dtype=float)
    ok = np.isfinite(N) & np.isfinite(X)
    N, X = N[ok], X[ok]
    base = dict(variable=variable, form=form, N_target=float(N_target),
                X_now=float(X_now))

    if N.size < config.min_points:
        return PhaseExtrapolation(accepted=False, reasons=[
            f"only {N.size} points in window (< {config.min_points})"], **base)
    if np.ptp(N) <= 0 or np.ptp(X) <= 0:
        return PhaseExtrapolation(accepted=False, reasons=[
            "no variation in N or X over the window"], **base)

    reasons: List[str] = []
    ratio = _extrapolation_ratio(N, N_target)
    if ratio > config.max_extrapolation_ratio:
        reasons.append(f"target N={N_target:+.2f} lies {ratio:.1f}x the sampled "
                       f"N-range beyond the window (> {config.max_extrapolation_ratio:g})")

    if form == "hyperbolic":
        if np.any(X <= 0):
            return PhaseExtrapolation(accepted=False, reasons=[
                "non-positive values: hyperbolic form needs X > 0"], **base)
        fit = fit_hyperbolic(N, X)
        fits = {"hyperbolic": fit}
        preds = {"hyperbolic": float(fit.predict(N_target))}
        a, b = fit.params
        if abs(fit.corr) < config.min_abs_corr:
            reasons.append(f"|corr(N, 1/{variable})|={abs(fit.corr):.2f} < "
                           f"{config.min_abs_corr}: not conduction-limited growth")
        if not np.isfinite(preds["hyperbolic"]):
            reasons.append(f"target N={N_target:+.2f} is at or past the asymptote "
                           f"a={a:+.2f}: {variable} would be infinite")
    else:
        lin, sat = fit_linear(N, X), fit_saturating(N, X)
        fits = {"linear": lin, "saturating": sat}
        preds = {"linear": float(lin.predict(N_target)),
                 "saturating": float(sat.predict(N_target))}
        if abs(lin.corr) < config.min_abs_corr:
            reasons.append(f"|corr({variable}, N)|={abs(lin.corr):.2f} < "
                           f"{config.min_abs_corr}: not slaved to the imbalance")
        scale = max(abs(preds["linear"] - X_now), abs(preds["saturating"] - X_now))
        if scale > 0:
            rel = abs(preds["linear"] - preds["saturating"]) / scale
            if rel > config.max_relative_disagreement:
                reasons.append(f"linear and saturating endpoints disagree by "
                               f"{rel:.0%} of the step: possible kink / "
                               f"feedback threshold")

    return PhaseExtrapolation(
        accepted=not reasons, reasons=reasons,
        prediction=preds[form] if not reasons else float("nan"),
        fits=fits, predictions=preds, extrapolation_ratio=ratio, **base)
