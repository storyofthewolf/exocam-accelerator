import numpy as np
import pytest

from exocam_accelerate.hindcast import (
    GateCalibration,
    annual_mean_series,
    calibrate_gate,
    default_origins,
    error_by_dt,
    hindcast_at,
    max_safe_dt,
    sweep_hindcasts,
)
from exocam_accelerate.safeguards import GateConfig
from exocam_accelerate.trends import TrendSeries


# --------------------------------------------------------------------------
# Series builders mimicking exocam-trend monthly output
# --------------------------------------------------------------------------

def monthly_linear(n_years=30, slope=0.5, intercept=250.0, seasonal=0.0, noise=0.0,
                   seed=0):
    """A monthly series X(t) = intercept + slope*t (+ seasonal + noise)."""
    rng = np.random.default_rng(seed)
    months = np.arange(1, n_years * 12 + 1)
    t = months / 12.0
    v = intercept + slope * t
    if seasonal:
        v = v + seasonal * np.sin(2 * np.pi * t)
    if noise:
        v = v + rng.normal(0, noise, size=t.size)
    return TrendSeries(times=t, values=v)


def monthly_exp_spinup(n_years=150, Teq=250.0, T0=230.0, tau=40.0, seed=1):
    """Monthly exponential relaxation toward equilibrium (the spin-up shape)."""
    months = np.arange(1, n_years * 12 + 1)
    t = months / 12.0
    v = Teq - (Teq - T0) * np.exp(-t / tau)
    return TrendSeries(times=t, values=v)


# --------------------------------------------------------------------------
# annual_mean_series
# --------------------------------------------------------------------------

class TestAnnualMean:
    def test_collapses_months_to_years(self):
        s = monthly_linear(n_years=10, slope=0.0, intercept=250.0, seasonal=5.0)
        ann = annual_mean_series(s)
        assert ann.values.shape == (10,)
        # seasonal cycle averages out over a full year
        assert np.allclose(ann.values, 250.0, atol=1e-6)

    def test_midyear_stamping(self):
        ann = annual_mean_series(monthly_linear(n_years=5))
        assert np.allclose(ann.times, [0.5, 1.5, 2.5, 3.5, 4.5])

    def test_drops_partial_trailing_year(self):
        # 30 months = 2 full years + 6 months; the partial year is dropped
        months = np.arange(1, 31)
        s = TrendSeries(times=months / 12.0, values=np.ones(30))
        ann = annual_mean_series(s)
        assert ann.values.shape == (2,)

    def test_preserves_linear_slope(self):
        s = monthly_linear(n_years=20, slope=0.5, intercept=100.0)
        ann = annual_mean_series(s)
        # annual means of a line lie on the same line at the mid-year points
        fit = np.polyfit(ann.times, ann.values, 1)
        assert np.isclose(fit[0], 0.5, atol=1e-6)

    def test_per_layer_series(self):
        months = np.arange(1, 25)
        t = months / 12.0
        vals = np.stack([250.0 + 0.1 * t, 240.0 - 0.2 * t], axis=1)  # (24, 2)
        ann = annual_mean_series(TrendSeries(times=t, values=vals))
        assert ann.values.shape == (2, 2)

    def test_too_few_full_years_raises(self):
        months = np.arange(1, 13)  # only one full year
        s = TrendSeries(times=months / 12.0, values=np.ones(12))
        with pytest.raises(ValueError):
            annual_mean_series(s)


# --------------------------------------------------------------------------
# hindcast_at
# --------------------------------------------------------------------------

class TestHindcastAt:
    def test_perfect_on_pure_line(self):
        # a noiseless line: forward-Euler extrapolation is exact
        ann = annual_mean_series(monthly_linear(n_years=40, slope=0.5))
        r = hindcast_at(ann, origin_year=20.5, window_years=10.0,
                        dt_years=10.0, variable="TS", max_abs_step=50.0)
        assert r is not None and r.accepted
        assert abs(r.error) < 1e-6
        # persistence would have been off by slope*dt = 5.0 (naive_error is
        # the max-abs magnitude, always non-negative)
        assert np.isclose(r.naive_error, 5.0, atol=1e-6)
        assert r.skill_score > 0.999

    def test_infeasible_target_past_end_returns_none(self):
        ann = annual_mean_series(monthly_linear(n_years=30, slope=0.5))
        r = hindcast_at(ann, origin_year=25.5, window_years=10.0,
                        dt_years=50.0, variable="TS", max_abs_step=50.0)
        assert r is None

    def test_window_uses_available_samples_near_start(self):
        # origin=3.5 with a 10-yr nominal window still has 4 samples (years
        # 0-3) in-window on a series starting at 0.5; that meets the gate's
        # sample floor, so the hindcast runs rather than being skipped.
        ann = annual_mean_series(monthly_linear(n_years=30, slope=0.5))
        r = hindcast_at(ann, origin_year=3.5, window_years=10.0,
                        dt_years=5.0, variable="TS", max_abs_step=50.0)
        assert r is not None and r.n_fit_points == 4

    def test_infeasible_when_too_few_samples_before_origin(self):
        # origin=1.5 leaves only years 0,1 before it -> below 2*n_segments
        ann = annual_mean_series(monthly_linear(n_years=30, slope=0.5))
        r = hindcast_at(ann, origin_year=1.5, window_years=10.0,
                        dt_years=5.0, variable="TS", max_abs_step=50.0)
        assert r is None

    def test_short_window_for_gate_returns_none(self):
        # window too short for the gate's 2*n_segments sample floor -> skip
        ann = annual_mean_series(monthly_linear(n_years=30, slope=0.5))
        r = hindcast_at(ann, origin_year=15.5, window_years=2.0,
                        dt_years=5.0, variable="TS", max_abs_step=50.0,
                        gate_config=GateConfig(n_segments=2))
        assert r is None

    def test_refusal_populates_baseline_not_prediction(self):
        # strong curvature -> gate refuses; predicted None, actual/naive filled
        ann = annual_mean_series(monthly_exp_spinup(n_years=60, tau=15.0))
        r = hindcast_at(ann, origin_year=8.5, window_years=6.0,
                        dt_years=20.0, variable="TS", max_abs_step=50.0,
                        gate_config=GateConfig(max_curvature_ratio=0.05,
                                               allow_decelerating=False))
        assert r is not None
        assert not r.accepted
        assert r.predicted is None and np.isnan(r.error)
        assert r.actual is not None and not np.isnan(r.naive_error)
        assert r.gate_reasons

    def test_clip_flag_reported(self):
        # huge slope over huge dt hits the clip
        ann = annual_mean_series(monthly_linear(n_years=60, slope=2.0))
        r = hindcast_at(ann, origin_year=30.5, window_years=10.0,
                        dt_years=25.0, variable="TS", max_abs_step=5.0)
        assert r is not None and r.accepted
        assert r.clipped

    def test_exp_spinup_error_has_expected_sign(self):
        # forward-Euler over a decaying exponential overshoots equilibrium:
        # predicted rises past actual, so error > 0 on a warming spin-up
        ann = annual_mean_series(monthly_exp_spinup())
        r = hindcast_at(ann, origin_year=20.5, window_years=10.0,
                        dt_years=40.0, variable="TS", max_abs_step=50.0)
        assert r is not None and r.accepted
        assert r.error > 0


# --------------------------------------------------------------------------
# default_origins
# --------------------------------------------------------------------------

class TestDefaultOrigins:
    def test_spacing_and_window_room(self):
        ann = annual_mean_series(monthly_linear(n_years=100))
        origins = default_origins(ann, window_years=10.0, step_years=20.0)
        # first origin leaves room for a full trailing window
        assert origins[0] >= ann.times[0] + 10.0 - 1.0
        assert np.allclose(np.diff(origins), 20.0)
        assert origins[-1] <= np.floor(ann.times[-1]) + 1e-9

    def test_empty_when_window_exceeds_span(self):
        ann = annual_mean_series(monthly_linear(n_years=5))
        assert default_origins(ann, window_years=20.0) == []


# --------------------------------------------------------------------------
# sweep_hindcasts
# --------------------------------------------------------------------------

class TestSweep:
    def test_sweep_skips_infeasible_and_keeps_feasible(self):
        s = monthly_linear(n_years=60, slope=0.3)
        res = sweep_hindcasts(
            s, "TS", max_abs_step=50.0,
            origins=[20.5, 40.5], windows=[10.0], dts=[5.0, 200.0],
        )
        # dt=200 is always infeasible on 60 yr; dt=5 feasible at both origins
        assert all(r.dt_years == 5.0 for r in res)
        assert len(res) == 2

    def test_already_annual_flag(self):
        ann = annual_mean_series(monthly_linear(n_years=60, slope=0.3))
        res = sweep_hindcasts(
            ann, "TS", max_abs_step=50.0,
            origins=[30.5], windows=[10.0], dts=[10.0],
            already_annual=True,
        )
        assert len(res) == 1 and res[0].accepted

    def test_pure_line_sweep_is_all_accurate(self):
        s = monthly_linear(n_years=80, slope=0.4)
        res = sweep_hindcasts(
            s, "TS", max_abs_step=50.0,
            origins=default_origins(annual_mean_series(s), 10.0, 10.0),
            windows=[10.0], dts=[10.0, 20.0],
        )
        assert res
        assert all(abs(r.error) < 1e-6 for r in res if r.accepted)


# --------------------------------------------------------------------------
# Aggregations: gate calibration, safe-dt, skill maps
# --------------------------------------------------------------------------

class TestAggregation:
    def test_calibrate_gate_counts(self):
        # hand-built results: 2 accepted-good, 1 accepted-bad, 1 refused-good,
        # 1 refused-bad
        def mk(accepted, error, naive):
            return type("R", (), {
                "accepted": accepted, "error": error, "naive_error": naive,
            })()
        results = [
            mk(True, 0.5, 3.0),    # accepted, |err|<=1 -> good
            mk(True, 0.9, 3.0),    # accepted good
            mk(True, 4.0, 3.0),    # accepted bad
            mk(False, np.nan, 0.5),  # refused, persistence within tol -> good (false alarm)
            mk(False, np.nan, 9.0),  # refused, persistence bad -> correct save
        ]
        cal = calibrate_gate(results, abs_tol=1.0)
        assert (cal.accepted_good, cal.accepted_bad) == (2, 1)
        assert (cal.refused_good, cal.refused_bad) == (1, 1)
        assert np.isclose(cal.miss_rate, 1 / 3)
        assert np.isclose(cal.false_alarm_rate, 1 / 2)

    def test_calibration_rates_nan_when_empty(self):
        cal = GateCalibration(tol=1.0, accepted_good=0, accepted_bad=0,
                              refused_good=0, refused_bad=0)
        assert np.isnan(cal.miss_rate) and np.isnan(cal.false_alarm_rate)

    def test_max_safe_dt_monotone_picks_largest_within_tol(self):
        # error grows with dt on the exp spin-up; there is a largest safe dt
        s = monthly_exp_spinup()
        res = sweep_hindcasts(
            s, "TS", max_abs_step=100.0,
            origins=default_origins(annual_mean_series(s), 10.0, 10.0),
            windows=[10.0], dts=[10.0, 25.0, 50.0],
        )
        safe = max_safe_dt(res, abs_tol=3.0)
        # small dt safe, largest dt not: the pick is a real dt from the set
        assert safe in (10.0, 25.0)
        # tightening tolerance can only shrink (or hold) the safe dt
        tighter = max_safe_dt(res, abs_tol=1.0)
        assert np.isnan(tighter) or tighter <= safe

    def test_max_safe_dt_nan_when_nothing_within_tol(self):
        s = monthly_exp_spinup()
        res = sweep_hindcasts(
            s, "TS", max_abs_step=100.0,
            origins=default_origins(annual_mean_series(s), 10.0, 10.0),
            windows=[10.0], dts=[50.0],
        )
        assert np.isnan(max_safe_dt(res, abs_tol=1e-6))

    def test_error_by_dt_structure_and_monotonicity(self):
        s = monthly_exp_spinup()
        res = sweep_hindcasts(
            s, "TS", max_abs_step=100.0,
            origins=default_origins(annual_mean_series(s), 10.0, 10.0),
            windows=[10.0], dts=[10.0, 25.0, 50.0],
        )
        table = error_by_dt(res)
        assert set(table) == {10.0, 25.0, 50.0}
        for stats in table.values():
            assert {"n", "mean", "median", "p90", "max"} <= set(stats)
        # mean error grows with dt for exponential spin-up
        means = [table[dt]["mean"] for dt in (10.0, 25.0, 50.0)]
        assert means[0] < means[1] < means[2]


# --------------------------------------------------------------------------
# Multi-element (per-layer / pointwise) series flow through unchanged
# --------------------------------------------------------------------------

class TestMultiElement:
    def test_per_layer_hindcast_reduces_to_maxabs(self):
        months = np.arange(1, 40 * 12 + 1)
        t = months / 12.0
        # two layers with different linear slopes
        vals = np.stack([250.0 + 0.5 * t, 240.0 + 0.1 * t], axis=1)
        s = TrendSeries(times=t, values=vals)
        ann = annual_mean_series(s)
        r = hindcast_at(ann, origin_year=20.5, window_years=10.0,
                        dt_years=10.0, variable="Tlev", max_abs_step=50.0)
        assert r is not None and r.accepted
        assert r.error_full.shape == (2,)
        # exact on lines -> full error ~0 in every layer
        assert np.allclose(r.error_full, 0.0, atol=1e-6)
