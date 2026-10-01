"""Pre-step safeguards: the trustworthiness gate and the hard magnitude clip.

The gate embodies the central safety principle of this tool: extrapolation
is only valid while the system is drifting toward equilibrium. Curvature
(d²X/dt²) that *reinforces* the trend (an accelerating drift), or a tendency
whose sign significantly reverses, signals proximity to a nonlinear
climate-feedback threshold (ice-albedo instability, runaway greenhouse) —
there the correct action is to REFUSE the step outright, not to shrink it
blindly: a smaller step across a bifurcation is still a step across a
bifurcation. Curvature that *opposes* the trend (a decelerating, asymptotic
approach) is the normal shape of a converging run and passes; sign flips
inside the noise are noise (user decision 2026-10-01; ``GateConfig``
``allow_decelerating`` / ``sign_flip_sigmas`` restore the original rules).

The clip is the independent second line of defense: an absolute per-step
magnitude limit no extrapolation may exceed, whatever the trend says
(Turbet et al. 2021 clip ΔT at 50 K per step for numerical stability).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .trends import TrendSeries, fit_curvature, fit_tendency, segment_slopes


@dataclass(frozen=True)
class GateConfig:
    """Tunables for the trustworthiness gate.

    max_curvature_ratio: refuse when |curvature contribution| over the window
        exceeds this fraction of the |linear change| over the window. With
        x(t) = a + bt + ct², the ratio is |c·T²| / max(|b·T|, negligible_change).
    negligible_change: linear change (variable units) below which a layer is
        treated as converged — it passes the gate with ~zero tendency rather
        than tripping the curvature ratio's denominator. Default 0.0 is the
        most conservative setting: any curvature atop a flat trend refuses.
    n_segments: sub-windows for the sign-stability check; the fitted slope in
        every segment must agree in sign with the full-window slope (segments
        whose |change| is below negligible_change are treated as neutral).
    """

    max_curvature_ratio: float = 0.5
    negligible_change: float = 0.0
    n_segments: int = 2
    #: curvature that *opposes* the trend (|dX/dt| shrinking over the window:
    #: a decelerating, asymptotic approach to equilibrium) is the normal shape
    #: of a converging run, not a threshold signal; only curvature that
    #: *reinforces* the trend (accelerating) is refused. False restores the
    #: original rule (any large curvature refuses).
    allow_decelerating: bool = True
    #: a sub-window counts against sign stability only if its slope disagrees
    #: with the full-window slope by more than this many standard errors
    #: (interannual noise flips short-segment slopes on a noisy, slow run).
    #: 0 restores the original rule (any disagreeing sign refuses).
    sign_flip_sigmas: float = 2.0


@dataclass(frozen=True)
class GateResult:
    """Per-element verdicts plus the overall decision.

    The gate is conservative: ``ok`` is True only if every element (layer,
    or gridpoint for pointwise use) passes both checks. Callers wanting
    per-layer selectivity can inspect ``ok_mask``, but the default stepper
    refuses the whole variable on any failure.
    """

    ok: bool
    ok_mask: np.ndarray
    curvature_ratio: np.ndarray
    sign_stable: np.ndarray
    converged_mask: np.ndarray
    reasons: tuple = field(default=())

    def __bool__(self) -> bool:
        return self.ok


def assess_trustworthiness(series: TrendSeries, config: GateConfig = GateConfig()) -> GateResult:
    """Decide whether a trend window supports linear extrapolation.

    Two checks, both required (Wordsworth/Turbet generalization):

    1. Curvature: the quadratic term's contribution over the window must be
       small relative to the linear change (see GateConfig.max_curvature_ratio).
    2. Sign stability: the tendency must not change sign across sub-windows.

    Elements whose total linear change is below ``negligible_change`` count as
    converged and pass automatically (their extrapolation step is ~zero).
    """
    T = series.window_length
    b, d2 = fit_curvature(series)
    c = 0.5 * d2

    linear_change = np.abs(b) * T
    quad_change = np.abs(c) * T * T
    converged = linear_change <= config.negligible_change

    denom = np.maximum(linear_change, config.negligible_change)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(denom > 0, quad_change / np.where(denom > 0, denom, 1.0), np.inf)
    ratio = np.where(quad_change == 0, 0.0, ratio)
    # b is the slope at the window midpoint, c the quadratic coefficient:
    # opposite signs = |slope| shrinking over the window (decelerating)
    decelerating = (np.sign(c) * np.sign(b)) < 0
    big_curv = ~(converged | (ratio <= config.max_curvature_ratio))
    curvature_ok = ~big_curv | (decelerating if config.allow_decelerating else False)

    slopes, slope_se = _segment_slopes_se(series, config.n_segments)
    full_slope = fit_tendency(series)
    seg_lengths = np.diff(np.linspace(series.times[0], series.times[-1],
                                      config.n_segments + 1))
    seg_change = np.abs(slopes) * seg_lengths.reshape(-1, *([1] * full_slope.ndim))
    seg_neutral = seg_change <= config.negligible_change
    # a disagreeing segment counts only if its slope is significantly of the
    # other sign (its own standard error)
    seg_noise = np.abs(slopes) <= config.sign_flip_sigmas * slope_se
    seg_agrees = (np.sign(slopes) == np.sign(full_slope)) | seg_neutral | seg_noise
    sign_stable = converged | seg_agrees.all(axis=0)

    ok_mask = curvature_ok & sign_stable
    ok = bool(np.all(ok_mask))

    reasons = []
    if not np.all(curvature_ok):
        n = int(np.sum(~curvature_ok))
        reasons.append(
            f"accelerating curvature: ratio exceeds {config.max_curvature_ratio} in "
            f"{n} element(s) (max {float(np.nanmax(ratio)):.3g}) and the trend is "
            f"speeding up — possible nonlinear feedback threshold; refusing to "
            f"extrapolate"
        )
    if not np.all(sign_stable):
        n = int(np.sum(~sign_stable))
        reasons.append(
            f"tendency sign significantly reversed across {config.n_segments} "
            f"sub-windows in {n} element(s); refusing to extrapolate"
        )

    return GateResult(
        ok=ok,
        ok_mask=np.asarray(ok_mask),
        curvature_ratio=np.asarray(ratio),
        sign_stable=np.asarray(sign_stable),
        converged_mask=np.asarray(converged),
        reasons=tuple(reasons),
    )


def _segment_slopes_se(series: TrendSeries, n_segments: int):
    """Per-segment linear slopes (as ``segment_slopes``) and their standard
    errors from each segment's residual scatter; shape (n_segments, *state)."""
    slopes = segment_slopes(series, n_segments)
    nt = len(series.times)
    bounds = np.linspace(0, nt, n_segments + 1).astype(int)
    se = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        t = series.times[lo:hi]
        y = series.values[lo:hi].reshape(hi - lo, -1)
        n = hi - lo
        if n < 3:
            se.append(np.zeros(y.shape[1]))
            continue
        tc = t - t.mean()
        sxx = float((tc ** 2).sum())
        A = np.vstack([np.ones(n), tc]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ coef
        se.append(np.sqrt((res ** 2).sum(0) / (n - 2) / sxx))
    return slopes, np.stack(se).reshape(slopes.shape)


@dataclass(frozen=True)
class ClipReport:
    """What the hard clip actually did."""

    clipped: np.ndarray
    n_clipped: int
    max_requested: float
    limit: float

    @property
    def any_clipped(self) -> bool:
        return self.n_clipped > 0


def clip_step(delta, max_abs_step: float) -> ClipReport:
    """Hard per-step magnitude clip, elementwise and sign-preserving.

    ``max_abs_step`` is an absolute limit in the variable's units,
    configurable per variable (Turbet: 50 K for temperature). It applies
    regardless of what the extrapolation suggested.
    """
    if max_abs_step <= 0:
        raise ValueError("max_abs_step must be positive")
    delta = np.asarray(delta, dtype=float)
    clipped = np.clip(delta, -max_abs_step, max_abs_step)
    n = int(np.sum(np.abs(delta) > max_abs_step))
    max_req = float(np.nanmax(np.abs(delta))) if delta.size else 0.0
    return ClipReport(clipped=clipped, n_clipped=n,
                      max_requested=max_req, limit=float(max_abs_step))
