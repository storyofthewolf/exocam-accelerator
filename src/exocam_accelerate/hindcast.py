"""Tier-0 offline hindcast validation.

The cheapest and first validation in the test plan (docs §6b, Tier 0): it needs
no new model runs and touches no restart files. Given an archived time series
that already ran to (near) equilibrium, we ask the question the whole tool
rests on — *would a forward-Euler tendency extrapolation have predicted where
the simulation actually went?* — and score the answer against ground truth the
archive already contains.

The unit operation is a single hindcast:

  1. Pick an origin time N inside the series.
  2. Fit a tendency over the trailing window [N-W, N] (via the real
     ``propose_step`` pipeline, gate and clip included).
  3. Extrapolate X(N) + delta out to N+dt.
  4. Compare to the *actual* archived value X(N+dt).

Sweeping (N, W, dt) over the 15 cold cases yields the three Tier-0
deliverables:

  * skill maps — forward-Euler error per variable / regime / dt;
  * gate calibration — of the (N, dt) pairs the gate *accepted*, how many the
    extrapolation actually got right, and how many it *refused* that would in
    fact have been fine (false-alarm rate). The gate earns its keep by
    refusing exactly the pairs a hindcast would have missed;
  * a defensible dt schedule — the largest dt whose accepted hindcast error
    stays within tolerance, as a function of how far the case is from
    equilibrium.

This module is pure computation over TrendSeries. It reads exocam-trend text
only through ``trend_io`` (the sole file interface); it imports no netCDF and
writes no restart files. It runs identically on the HPC or a laptop — the
input is a few MB of text.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

import numpy as np

from .safeguards import GateConfig
from .stepper import propose_step
from .trends import TrendSeries


# --------------------------------------------------------------------------
# Reducing a raw monthly series to the annual-mean series the hindcast scores
# --------------------------------------------------------------------------

def annual_mean_series(series: TrendSeries, months_per_year: int = 12) -> TrendSeries:
    """Collapse a monthly TrendSeries to one value per whole year.

    docs §6a(4): tendencies must be measured from annual means, never raw
    monthly values, or Δt multiplies the seasonal cycle. Each calendar year is
    averaged and stamped at its mid-year (year + 0.5); any partial year at
    either end (fewer than ``months_per_year`` samples) is dropped.

    exocam-trend numbers months from 1, writing ``month`` as the running index
    (so December of the first year is month 12 -> time 1.0 year). Binning must
    therefore group months (1..12) -> year 0, (13..24) -> year 1, i.e. by
    ``(month_index - 1) // 12`` — not by ``floor(time)``, which would split
    each December off into the next year's bin.

    Works for scalar or per-layer/pointwise series (any trailing state shape).
    """
    if series.times.size == 0:
        raise ValueError("empty series")
    # recover the 1-based running month index from the year-valued times
    month_idx = np.rint(series.times * months_per_year).astype(int)
    year = (month_idx - 1) // months_per_year

    uniq, counts = np.unique(year, return_counts=True)
    full = counts == months_per_year
    keep = uniq[full]
    if len(keep) < 2:
        raise ValueError(
            "need at least 2 complete years after annual averaging "
            f"(found {len(keep)} full year(s) of {months_per_year} samples)"
        )
    times = keep.astype(float) + 0.5
    vals = np.stack([series.values[year == y].mean(axis=0) for y in keep])
    return TrendSeries(times=times, values=vals)


# --------------------------------------------------------------------------
# One hindcast
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class HindcastResult:
    """Outcome of one (origin, window, dt) hindcast for one variable.

    ``accepted`` is the gate verdict. When accepted, ``predicted`` is
    X(origin) + clipped delta and ``error = predicted - actual``; when refused,
    those extrapolation fields are None but ``actual``/``naive_error`` are still
    filled so a refusal can be scored against what persistence would have done.
    """

    variable: str
    origin_year: float
    window_years: float
    dt_years: float
    n_fit_points: int

    accepted: bool
    gate_reasons: tuple

    actual: Optional[np.ndarray]        # X(origin + dt), ground truth
    baseline: Optional[np.ndarray]      # X(origin), the "do nothing" prediction
    predicted: Optional[np.ndarray]     # X(origin) + delta (None if refused)

    # scalar summaries (np.nan where undefined); for multi-element series these
    # are the max-abs over elements, with full arrays in the *_full fields.
    error: float                        # predicted - actual (accepted only)
    naive_error: float                  # baseline - actual (always defined)
    skill_score: float                  # 1 - |error| / |naive_error|

    error_full: Optional[np.ndarray] = None
    naive_error_full: Optional[np.ndarray] = None
    clipped: bool = False


def _reduce(x: Optional[np.ndarray]) -> float:
    if x is None:
        return float("nan")
    a = np.abs(np.asarray(x, dtype=float))
    return float(np.nanmax(a)) if a.size else float("nan")


def hindcast_at(
    series: TrendSeries,
    origin_year: float,
    window_years: float,
    dt_years: float,
    variable: str,
    max_abs_step: float,
    gate_config: GateConfig = GateConfig(),
    *,
    time_tol: float = 0.51,
) -> Optional[HindcastResult]:
    """Run one hindcast on an already annual-mean series, or None if infeasible.

    Returns None (rather than raising) when the series cannot support this
    (origin, window, dt) triple — no fit window, or no archived sample near
    origin+dt to score against — so a sweep can skip cleanly.

    ``time_tol`` (years) is how close an archived sample must sit to the
    requested origin / target time to count as a match; the default half-year
    tolerance matches the mid-year stamping of ``annual_mean_series``.
    """
    t = series.times

    # locate the fit window [origin - window, origin]
    lo, hi = origin_year - window_years, origin_year
    in_win = (t >= lo - time_tol) & (t <= hi + time_tol)
    # The gate's sign-stability check needs >= 2*n_segments samples; a window
    # with fewer cannot be assessed, so the triple is infeasible (skip, don't
    # crash). Curvature needs >= 3; the segment bound is the stricter one.
    if int(np.sum(in_win)) < 2 * gate_config.n_segments:
        return None
    win = TrendSeries(times=t[in_win], values=series.values[in_win])

    # baseline (value at origin) and target (value at origin + dt)
    i_origin = int(np.argmin(np.abs(t - origin_year)))
    if abs(t[i_origin] - origin_year) > time_tol:
        return None
    target_time = origin_year + dt_years
    i_target = int(np.argmin(np.abs(t - target_time)))
    if abs(t[i_target] - target_time) > time_tol:
        return None

    baseline = series.values[i_origin]
    actual = series.values[i_target]
    naive_full = baseline - actual

    proposal = propose_step(win, dt_years, max_abs_step, gate_config)

    if not proposal.accepted:
        return HindcastResult(
            variable=variable,
            origin_year=float(origin_year),
            window_years=float(window_years),
            dt_years=float(dt_years),
            n_fit_points=len(win.times),
            accepted=False,
            gate_reasons=proposal.gate.reasons,
            actual=np.asarray(actual),
            baseline=np.asarray(baseline),
            predicted=None,
            error=float("nan"),
            naive_error=_reduce(naive_full),
            skill_score=float("nan"),
            error_full=None,
            naive_error_full=np.asarray(naive_full),
        )

    predicted = baseline + proposal.delta
    err_full = predicted - actual
    err = _reduce(err_full)
    naive = _reduce(naive_full)
    skill = 1.0 - (err / naive) if naive > 0 else float("nan")

    return HindcastResult(
        variable=variable,
        origin_year=float(origin_year),
        window_years=float(window_years),
        dt_years=float(dt_years),
        n_fit_points=len(win.times),
        accepted=True,
        gate_reasons=(),
        actual=np.asarray(actual),
        baseline=np.asarray(baseline),
        predicted=np.asarray(predicted),
        error=err,
        naive_error=naive,
        skill_score=skill,
        error_full=np.asarray(err_full),
        naive_error_full=np.asarray(naive_full),
        clipped=bool(proposal.clip.any_clipped),
    )


# --------------------------------------------------------------------------
# The sweep
# --------------------------------------------------------------------------

def sweep_hindcasts(
    series: TrendSeries,
    variable: str,
    max_abs_step: float,
    origins: Sequence[float],
    windows: Sequence[float],
    dts: Sequence[float],
    gate_config: GateConfig = GateConfig(),
    *,
    already_annual: bool = False,
    time_tol: float = 0.51,
) -> List[HindcastResult]:
    """Every feasible (origin, window, dt) hindcast for one variable/series.

    Unless ``already_annual``, the input monthly series is reduced with
    ``annual_mean_series`` first. Infeasible triples (window off the start,
    target past the end) are silently skipped.
    """
    ann = series if already_annual else annual_mean_series(series)
    out: List[HindcastResult] = []
    for w in windows:
        for n in origins:
            for dt in dts:
                r = hindcast_at(
                    ann, n, w, dt, variable, max_abs_step,
                    gate_config, time_tol=time_tol,
                )
                if r is not None:
                    out.append(r)
    return out


def default_origins(series: TrendSeries, window_years: float,
                     step_years: float = 5.0) -> List[float]:
    """Origin years spanning the archive with room for the fit window.

    Every ``step_years`` from the earliest year that admits a full trailing
    window up to the last archived year. On an annual-mean series pass the
    annualised series; on a raw monthly one the year span is the same.
    """
    y0 = float(np.floor(series.times[0]))
    y1 = float(np.floor(series.times[-1]))
    first = y0 + window_years
    if first > y1:
        return []
    n = int(np.floor((y1 - first) / step_years))
    return [first + k * step_years for k in range(n + 1)]


# --------------------------------------------------------------------------
# Aggregation into the three Tier-0 deliverables
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class GateCalibration:
    """How well the gate's accept/refuse decision matched hindcast outcomes.

    A hindcast is "good" when its accepted extrapolation beat persistence by
    the tolerance (skill within band), "bad" otherwise. Against that ground
    truth:

      accepted_good   gate said go, extrapolation was right   (true accept)
      accepted_bad    gate said go, extrapolation missed       (MISS — costly)
      refused_good    gate refused, extrapolation would've been right
                                                               (false alarm)
      refused_bad     gate refused, extrapolation would've missed (correct save)

    The gate is well-calibrated when accepted_bad and refused_good are both
    small: it lets through what works and stops what doesn't. ``miss_rate`` and
    ``false_alarm_rate`` are the two numbers to drive threshold tuning.
    """

    tol: float
    accepted_good: int
    accepted_bad: int
    refused_good: int
    refused_bad: int

    @property
    def n_accepted(self) -> int:
        return self.accepted_good + self.accepted_bad

    @property
    def n_refused(self) -> int:
        return self.refused_good + self.refused_bad

    @property
    def miss_rate(self) -> float:
        """Fraction of accepted steps that missed — the dangerous quantity."""
        return self.accepted_bad / self.n_accepted if self.n_accepted else float("nan")

    @property
    def false_alarm_rate(self) -> float:
        """Fraction of refusals that would in fact have been fine."""
        return self.refused_good / self.n_refused if self.n_refused else float("nan")


def calibrate_gate(results: Iterable[HindcastResult], abs_tol: float) -> GateCalibration:
    """Score gate decisions against hindcast outcomes at an absolute tolerance.

    ``abs_tol`` is in the variable's own units (e.g. K for temperature, m for
    ice thickness): a hindcast is "good" if its absolute error would have been
    <= abs_tol. For refused steps we still know what the extrapolation *would*
    have produced only if we recompute it; here we instead use the honest
    fallback — a refusal is scored good/bad by whether *persistence* (doing
    nothing, the actual consequence of refusing) lands within tol. That makes
    refused_good = "we refused but even doing nothing was within tol, so the
    refusal cost us a safe acceleration opportunity."
    """
    ag = ab = rg = rb = 0
    for r in results:
        if r.accepted:
            if not np.isnan(r.error) and abs(r.error) <= abs_tol:
                ag += 1
            else:
                ab += 1
        else:
            # consequence of refusing is persistence; naive_error is its error
            if not np.isnan(r.naive_error) and abs(r.naive_error) <= abs_tol:
                rg += 1
            else:
                rb += 1
    return GateCalibration(tol=abs_tol, accepted_good=ag, accepted_bad=ab,
                           refused_good=rg, refused_bad=rb)


def max_safe_dt(results: Iterable[HindcastResult], abs_tol: float) -> float:
    """Largest dt whose *accepted* hindcasts all stayed within abs_tol.

    The evidence for a defensible decreasing-Δt schedule (docs §6b Tier-0
    deliverable c): scan accepted results by dt and return the largest dt at
    which no accepted extrapolation exceeded the tolerance. Returns nan if no
    accepted result meets the tolerance at any dt.
    """
    by_dt = {}
    for r in results:
        if not r.accepted or np.isnan(r.error):
            continue
        by_dt.setdefault(r.dt_years, []).append(abs(r.error))
    safe = [dt for dt, errs in by_dt.items() if max(errs) <= abs_tol]
    return max(safe) if safe else float("nan")


def error_by_dt(results: Iterable[HindcastResult]) -> dict:
    """Map dt -> summary of accepted-hindcast |error| (for skill maps).

    Each value is a dict with n, mean, median, p90, max of the absolute error
    over accepted hindcasts at that dt. The raw material for a skill map.
    """
    buckets = {}
    for r in results:
        if r.accepted and not np.isnan(r.error):
            buckets.setdefault(r.dt_years, []).append(abs(r.error))
    out = {}
    for dt, errs in sorted(buckets.items()):
        a = np.asarray(errs)
        out[dt] = {
            "n": int(a.size),
            "mean": float(a.mean()),
            "median": float(np.median(a)),
            "p90": float(np.percentile(a, 90)),
            "max": float(a.max()),
        }
    return out
