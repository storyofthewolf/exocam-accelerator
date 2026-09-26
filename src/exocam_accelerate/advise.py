"""Jump advisor: from a case's exocam-trend output to a recommended ice factor.

Pure computation over exocam-trend columns (read via ``trend_io``).

Physics (docs/phase-space-extrapolation.md, "Ice growth is Stefan-limited"):
in the cold cases the TOA deficit N is conducted through the sea ice and
freezes onto its base, so N*hi is nearly constant: N = a + b/hi. The ice never
"equilibrates" at N = 0 (hi diverges there); instead a jump chooses a target
imbalance and sets the ice to the thickness the conduction law says produces
it. That target is expressed as the fraction of the current imbalance to
remove (``n_fraction``), or an absolute ``N_target``.

The actionable output is the **ice factor** = hi_target / hi_now, hard-clipped,
which the restart layer applies to ``vicen``/``eicen`` in cice.r. Also
reported, as the reference for checking the post-jump run:

* ``N_after`` — the imbalance the conduction law predicts for the clipped
  factor; the post-jump run's N should settle near it within a few years.
* ``TS``/``Tsfc`` at ``N_after`` (linear in N, Gregory-style).
* ``years_skipped`` — the time the run would have needed to grow that ice,
  from the window's h^2-vs-t slope (Stefan: h^2 grows linearly in time).
* ``qi``/``hi`` drift — scaling vicen and eicen together assumes enthalpy per
  unit volume is steady.

Gate additions beyond ``phase_space``: the ice edge must have settled
(|dICEFRAC| over the window small). Scaling vicen with aicen held fixed, and
the Stefan law itself, both assume fixed ice area; windows still inside the
initial ice-edge transient were the only failures in the hindcast.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from .hindcast import annual_mean_series
from .phase_space import PhaseExtrapolation, PhaseGateConfig, extrapolate
from .trend_io import trend_series

ICE_VOLUME_VAR = "hi"
ICE_ENTHALPY_VAR = "qi"
ICE_AREA_VAR = "ICEFRAC"
TEMPERATURE_VARS = ("TS", "Tsfc")


@dataclass(frozen=True)
class AdvisorConfig:
    #: fit window in years; None = the longest of ``auto_windows`` over which
    #: the ice edge has settled (hindcast: accuracy is ~flat in window length,
    #: longer windows are accepted far more often late in a run)
    window_years: Optional[float] = None
    auto_windows: tuple = (40.0, 30.0, 20.0, 10.0)
    which: str = "int2"
    imbalance: str = "energy_top"
    #: fraction of the current imbalance the jump should remove
    n_fraction: float = 0.5
    #: absolute target imbalance (W/m2); overrides n_fraction when set
    N_target: Optional[float] = None
    max_ice_factor: float = 1.5
    max_icefrac_change: float = 0.02
    #: warn when qi/hi drifts by more than this fraction over the window
    max_enthalpy_drift: float = 0.05
    gate: PhaseGateConfig = field(default_factory=PhaseGateConfig)


@dataclass(frozen=True)
class Advice:
    case: str
    year_now: float           # mid-point of the latest complete model year
    N_now: float              # latest native annual mean
    N_now_fit: float          # conduction law at hi_now (noise-free)
    N_target: float
    config: AdvisorConfig
    window_years: float       # the fit window actually used
    ice: Optional[PhaseExtrapolation]
    temperatures: Dict[str, PhaseExtrapolation]
    ice_factor: Optional[float]          # None -> do not jump
    ice_factor_raw: Optional[float]
    ice_factor_clipped: bool
    N_after: Optional[float]
    years_skipped: Optional[float]
    reasons: tuple                        # why no jump (empty when jumping)
    warnings: tuple

    @property
    def jump(self) -> bool:
        return self.ice_factor is not None

    def to_dict(self) -> dict:
        def num(x):
            return None if x is None or not np.isfinite(x) else float(x)

        def ext(r: PhaseExtrapolation):
            return {
                "form": r.form,
                "accepted": r.accepted,
                "reasons": list(r.reasons),
                "N_target": num(r.N_target),
                "X_now": num(r.X_now),
                "prediction": num(r.prediction),
                "predictions": {k: num(v) for k, v in r.predictions.items()},
                "fits": {k: {"params": list(f.params), "corr": f.corr,
                             "rmse": f.rmse} for k, f in r.fits.items()},
                "extrapolation_ratio": num(r.extrapolation_ratio),
            }

        c = self.config
        return {
            "case": self.case,
            "year_now": self.year_now,
            "N_now": num(self.N_now),
            "N_now_fit": num(self.N_now_fit),
            "N_target": num(self.N_target),
            "config": {
                "window_years": self.window_years, "which": c.which,
                "imbalance": c.imbalance, "n_fraction": c.n_fraction,
                "N_target": c.N_target, "max_ice_factor": c.max_ice_factor,
                "max_icefrac_change": c.max_icefrac_change,
            },
            "ice_factor": self.ice_factor,
            "ice_factor_raw": num(self.ice_factor_raw),
            "ice_factor_clipped": self.ice_factor_clipped,
            "N_after": num(self.N_after),
            "years_skipped": num(self.years_skipped),
            "expected_after_jump": {v: num(r.prediction)
                                    for v, r in self.temperatures.items()},
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "ice": ext(self.ice) if self.ice is not None else None,
            "temperatures": {v: ext(r) for v, r in self.temperatures.items()},
        }


def _annual(columns, variable, which):
    s = annual_mean_series(trend_series(columns, variable, which=which))
    return s.times, s.values


def _have(columns, var, which) -> bool:
    return f"{var}_{which}" in columns and f"{var}_native" in columns


def _choose_window(columns, t, config: AdvisorConfig) -> float:
    """Fixed window, or the longest candidate with a settled ice edge."""
    if config.window_years is not None:
        return float(config.window_years)
    candidates = sorted(config.auto_windows, reverse=True)
    if not _have(columns, ICE_AREA_VAR, config.which):
        return float(candidates[-1])
    _, f = _annual(columns, ICE_AREA_VAR, config.which)
    for W in candidates:
        w = t > t[-1] - W
        if abs(f[w][-1] - f[w][0]) <= config.max_icefrac_change:
            return float(W)
    return float(candidates[-1])   # the gate below will refuse it


def advise(columns: Dict[str, np.ndarray], case: str,
           config: AdvisorConfig = AdvisorConfig()) -> Advice:
    """Build jump advice from merged exocam-trend columns (``load_case``)."""
    t, N_fit = _annual(columns, config.imbalance, config.which)
    _, N_nat = _annual(columns, config.imbalance, "native")
    year_now = float(t[-1])
    N_now = float(N_nat[-1])
    reasons, warnings = [], []
    window = _choose_window(columns, t, config)
    w = t > t[-1] - window

    def empty(N_target=float("nan"), N_now_fit=float("nan"), ice=None, temps=None):
        return Advice(case, year_now, N_now, N_now_fit, N_target, config, window,
                      ice, temps or {}, None, None, False, None, None,
                      tuple(reasons), tuple(warnings))

    def have(var):
        return _have(columns, var, config.which)

    if not have(ICE_VOLUME_VAR):
        reasons.append(f"no {ICE_VOLUME_VAR} series in the trend output")
        return empty()
    _, h_fit = _annual(columns, ICE_VOLUME_VAR, config.which)
    _, h_nat = _annual(columns, ICE_VOLUME_VAR, "native")
    h_now = float(h_nat[-1])

    # Fit once to locate the current point on the conduction law, then set the
    # target relative to it (the native N_now is interannually noisy).
    probe = extrapolate(ICE_VOLUME_VAR, N_fit[w], h_fit[w], N_now, h_now,
                        "hyperbolic", config.gate)
    if "hyperbolic" not in probe.fits:
        reasons.extend(probe.reasons)
        return empty(ice=probe)
    a, b = probe.fits["hyperbolic"].params
    N_now_fit = a + b / h_now
    if config.N_target is not None:
        N_target = float(config.N_target)
    else:
        N_target = N_now_fit * (1.0 - config.n_fraction)
    ice = extrapolate(ICE_VOLUME_VAR, N_fit[w], h_fit[w], N_target, h_now,
                      "hyperbolic", config.gate)
    reasons.extend(ice.reasons)

    if have(ICE_AREA_VAR):
        _, f_fit = _annual(columns, ICE_AREA_VAR, config.which)
        d_area = float(abs(f_fit[w][-1] - f_fit[w][0]))
        if d_area > config.max_icefrac_change:
            reasons.append(f"ice edge still moving: |d{ICE_AREA_VAR}|={d_area:.3f} "
                           f"over the window (> {config.max_icefrac_change}); the "
                           f"conduction law and fixed-area scaling need a settled edge")
    else:
        warnings.append(f"no {ICE_AREA_VAR} series: cannot confirm the ice edge settled")

    if have(ICE_ENTHALPY_VAR):
        _, q_fit = _annual(columns, ICE_ENTHALPY_VAR, config.which)
        ratio = q_fit[w] / h_fit[w]
        drift = float(abs(ratio[-1] - ratio[0]) / abs(ratio.mean()))
        if drift > config.max_enthalpy_drift:
            warnings.append(f"{ICE_ENTHALPY_VAR}/{ICE_VOLUME_VAR} drifted {drift:.1%} "
                            f"over the window: enthalpy per volume is not steady, so "
                            f"scaling eicen with vicen is approximate")

    factor = raw = N_after = years = None
    clipped = False
    if not reasons:
        raw = ice.prediction / h_now
        lo, up = 1.0 / config.max_ice_factor, config.max_ice_factor
        factor = float(min(max(raw, lo), up))
        clipped = factor != raw
        h_new = h_now * factor
        N_after = float(a + b / h_new)
        # Stefan: h^2 grows linearly in time; its slope converts ice to years.
        h2_slope = np.polyfit(t[w], h_fit[w] ** 2, 1)[0]
        if h2_slope > 0:
            years = float((h_new**2 - h_now**2) / h2_slope)

    temps = {}
    N_ref = N_after if N_after is not None else N_target
    for var in TEMPERATURE_VARS:
        if not have(var):
            continue
        _, X_fit = _annual(columns, var, config.which)
        _, X_nat = _annual(columns, var, "native")
        temps[var] = extrapolate(var, N_fit[w], X_fit[w], N_ref, float(X_nat[-1]),
                                 "linear", config.gate)
        if factor is not None and not temps[var].accepted:
            warnings.append(f"{var} reference not available: "
                            + "; ".join(temps[var].reasons))

    return Advice(case, year_now, N_now, float(N_now_fit), float(N_target), config,
                  window, ice, temps, factor, raw, clipped, N_after, years,
                  tuple(reasons), tuple(warnings))
