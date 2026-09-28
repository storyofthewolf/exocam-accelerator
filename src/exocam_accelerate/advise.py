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
reported, as the reference for checking the post-jump run (``check``):

* the fitted conduction law (a, b) and ``N_after`` — the imbalance the law
  predicts for the clipped factor;
* ``TS``/``Tsfc`` as linear functions of N (Gregory-style);
* ``years_skipped`` — the time the run would have needed to grow that ice,
  from the window's h^2-vs-t slope (Stefan: h^2 grows linearly in time).

Gate additions beyond ``phase_space``: the ice edge must have settled (small
ICEFRAC change across the window). Scaling vicen with aicen held fixed, and the
Stefan law itself, both assume fixed ice area.

After a jump (``since_year``, or detected automatically as a one-year step in
hi): only data from ``settle_years`` after the jump are used, and they are
fitted as native annual means — the int2 column is a 10-yr trailing mean and
straddles the jump for a decade. Native annual N is noisy, so the post-jump
gate asks for >= ``post_jump_min_years`` of data and a looser correlation
(hindcast: 20 native years at |corr| >= 0.7 match int2 accuracy, ~3 % median).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .hindcast import annual_mean_series
from .phase_space import PhaseExtrapolation, PhaseGateConfig, extrapolate
from .trend_io import trend_series

ICE_VOLUME_VAR = "hi"
ICE_ENTHALPY_VAR = "qi"
ICE_AREA_VAR = "ICEFRAC"
TEMPERATURE_VARS = ("TS", "Tsfc")

#: Version of the serialized ``Advice.to_dict()`` shape. Bump this whenever a
#: field is added, removed, or reinterpreted; ``runstate.preflight`` refuses
#: advice whose ``schema_version`` is missing or not in
#: ``runstate.KNOWN_ADVICE_SCHEMA_VERSIONS``.
ADVICE_SCHEMA_VERSION = "1"


@dataclass(frozen=True)
class AdvisorConfig:
    #: fit window in years; None = the longest of ``auto_windows`` over which
    #: the ice edge has settled (hindcast: accuracy is ~flat in window length,
    #: longer windows are accepted far more often late in a run)
    window_years: Optional[float] = None
    auto_windows: Tuple[float, ...] = (40.0, 30.0, 20.0, 10.0)
    which: str = "int2"
    imbalance: str = "energy_top"
    #: fraction of the current imbalance the jump should remove
    n_fraction: float = 0.5
    #: absolute target imbalance (W/m2); overrides n_fraction when set
    N_target: Optional[float] = None
    max_ice_factor: float = 1.5
    #: max |ICEFRAC change| across the window (3-yr edge means)
    max_icefrac_change: float = 0.02
    #: warn when qi/hi drifts by more than this fraction over the window
    max_enthalpy_drift: float = 0.05
    gate: PhaseGateConfig = field(default_factory=PhaseGateConfig)

    # --- after a previous jump ---
    #: first model year run from a jumped state; None = detect from hi
    since_year: Optional[int] = None
    detect_jumps: bool = True
    #: a one-year hi ratio above this marks a jump (natural growth is ~1-2 %/yr)
    jump_detect_ratio: float = 1.08
    settle_years: int = 2
    post_jump_min_years: int = 15
    post_jump_gate: PhaseGateConfig = field(
        default_factory=lambda: PhaseGateConfig(min_abs_corr=0.7))

    def __post_init__(self) -> None:
        if not (0.0 < self.n_fraction < 1.0):
            raise ValueError(f"n_fraction must satisfy 0 < n_fraction < 1, got "
                             f"{self.n_fraction!r}")
        if not (1.0 <= self.max_ice_factor <= 2.0):
            raise ValueError(f"max_ice_factor must satisfy 1 <= max_ice_factor <= 2 "
                             f"(a thickening-only cold-case jump), got "
                             f"{self.max_ice_factor!r}")
        if self.window_years is not None and self.window_years <= 0:
            raise ValueError(f"window_years must be positive, got "
                             f"{self.window_years!r}")
        if any(w <= 0 for w in self.auto_windows):
            raise ValueError(f"auto_windows entries must all be positive, got "
                             f"{self.auto_windows!r}")
        if self.settle_years <= 0:
            raise ValueError(f"settle_years must be positive, got "
                             f"{self.settle_years!r}")
        if self.post_jump_min_years <= 0:
            raise ValueError(f"post_jump_min_years must be positive, got "
                             f"{self.post_jump_min_years!r}")
        if self.max_icefrac_change < 0:
            raise ValueError("max_icefrac_change must be non-negative")
        if self.max_enthalpy_drift < 0:
            raise ValueError("max_enthalpy_drift must be non-negative")


@dataclass(frozen=True)
class Advice:
    case: str
    model_year: int           # latest complete model year in the trend data
    N_now: float              # its native annual mean
    N_now_fit: float          # conduction law at hi_now (noise-free)
    N_target: float
    config: AdvisorConfig
    which_used: str           # trend column actually fitted
    window_years: float       # the fit window actually used
    since_year: Optional[int]
    detected_jumps: Tuple[int, ...]
    ice: Optional[PhaseExtrapolation]
    temperatures: Dict[str, PhaseExtrapolation]
    ice_factor: Optional[float]          # None -> do not jump
    ice_factor_raw: Optional[float]
    ice_factor_clipped: bool
    N_after: Optional[float]
    years_skipped: Optional[float]
    reasons: tuple                        # why no jump (empty when jumping)
    warnings: tuple
    provenance: Optional[Dict[str, str]] = None  # input trend file -> sha256

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
            "schema_version": ADVICE_SCHEMA_VERSION,
            "case": self.case,
            "model_year": self.model_year,
            "N_now": num(self.N_now),
            "N_now_fit": num(self.N_now_fit),
            "N_target": num(self.N_target),
            "which": self.which_used,
            "window_years": self.window_years,
            "since_year": self.since_year,
            "detected_jumps": list(self.detected_jumps),
            "config": {
                "n_fraction": c.n_fraction, "N_target": c.N_target,
                "max_ice_factor": c.max_ice_factor, "imbalance": c.imbalance,
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
            "provenance": dict(self.provenance) if self.provenance else None,
        }


def _annual(columns, variable, which):
    s = annual_mean_series(trend_series(columns, variable, which=which))
    return s.times, s.values


def _have(columns, var, which) -> bool:
    return f"{var}_{which}" in columns and f"{var}_native" in columns


def model_years(t, start_year: int) -> np.ndarray:
    """Model year of each annual-mean bin (bins are stamped year-index + 0.5)."""
    return (start_year + np.floor(t)).astype(int)


def detect_jumps(years, hi, ratio: float = 1.08, contrast: float = 4.0,
                 neighbours: int = 5) -> List[int]:
    """Model years whose annual-mean hi steps up abruptly from the year before.

    A jump is a one-year growth ratio above ``ratio`` that also exceeds
    ``contrast`` x the median growth of the ``neighbours`` years on either
    side — early spin-up grows >8 %/yr for years on end, a jump stands alone.
    Returns the first model year run from each jumped state.
    """
    hi = np.asarray(hi, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = hi[1:] / hi[:-1] - 1.0
    out = []
    for i in np.where(g > ratio - 1.0)[0]:
        if i < neighbours:          # no history to contrast with (initial growth)
            continue
        around = np.r_[g[max(0, i - neighbours):i], g[i + 1:i + 1 + neighbours]]
        around = around[np.isfinite(around)]
        if around.size and g[i] > contrast * max(np.median(around), 0.0):
            out.append(int(years[i + 1]))
    return out


def _icefrac_change(t, f, mask, edge_years: int = 3) -> float:
    """|ICEFRAC change| across the window: mean of its last ``edge_years``
    minus mean of its first. An endpoint comparison (not a trend) so a window
    reaching back into the initial ice-edge transient is caught; the edge
    averaging keeps it robust on native annual data.
    """
    ff = f[mask]
    if ff.size < 2:
        return float("inf")
    k = max(1, min(edge_years, ff.size // 2))
    return float(abs(ff[-k:].mean() - ff[:k].mean()))


def advise(columns: Dict[str, np.ndarray], case: str,
           config: AdvisorConfig = AdvisorConfig(),
           start_year: int = 1,
           provenance: Optional[Dict[str, str]] = None) -> Advice:
    """Build jump advice from merged exocam-trend columns (``load_case``).

    ``start_year`` is the model year of the trend series' first month
    (``trend_io.case_start_year``); it only matters for reporting and for
    ``since_year``, which is a model year.

    ``provenance`` (e.g. ``trend_io.file_provenance(...)``) is recorded
    verbatim in the advice so it can be audited later; it is not otherwise
    used here.
    """
    t, N_nat = _annual(columns, config.imbalance, "native")
    years = model_years(t, start_year)
    model_year = int(years[-1])
    N_now = float(N_nat[-1])
    reasons: List[str] = []
    warnings: List[str] = []
    have_hi = _have(columns, ICE_VOLUME_VAR, "native")

    # ---- previous jumps: restrict to post-jump data, fit native ----
    detected: List[int] = []
    if have_hi and config.detect_jumps:
        _, h_nat_all = _annual(columns, ICE_VOLUME_VAR, "native")
        detected = detect_jumps(years, h_nat_all, config.jump_detect_ratio)
    since = config.since_year
    if since is None and detected:
        since = detected[-1]
        warnings.append(f"detected a jump in {ICE_VOLUME_VAR} at model year "
                        f"{since}; using post-jump data only (override with --since)")
    elif since is not None and detected and detected[-1] > since:
        warnings.append(f"a later jump is visible in {ICE_VOLUME_VAR} at model year "
                        f"{detected[-1]} than --since {since}")

    if since is not None:
        which, gate = "native", config.post_jump_gate
        usable = years >= since + config.settle_years
        n_usable = int(usable.sum())
        cap = float(n_usable)
    else:
        which, gate = config.which, config.gate
        usable = np.ones_like(t, dtype=bool)
        n_usable, cap = t.size, float("inf")

    _, N_fit = _annual(columns, config.imbalance, which)

    # ---- window ----
    area = None
    if _have(columns, ICE_AREA_VAR, which):
        _, area = _annual(columns, ICE_AREA_VAR, which)
    if config.window_years is not None:
        window = min(float(config.window_years), cap)
    else:
        cands = {W for W in config.auto_windows if W <= cap}
        if since is not None:
            # use all settled post-jump years (up to the longest candidate)
            cands.add(min(cap, max(config.auto_windows)))
            cands = {W for W in cands if W >= config.post_jump_min_years}
        cands = sorted(cands, reverse=True)
        window = cands[-1] if cands else cap
        if area is not None:
            for W in cands:
                m = usable & (t > t[-1] - W)
                if _icefrac_change(t, area, m) <= config.max_icefrac_change:
                    window = W
                    break
    w = usable & (t > t[-1] - window)

    def result(N_target=float("nan"), N_now_fit=float("nan"), ice=None, temps=None,
               factor=None, raw=None, clipped=False, N_after=None, skipped=None):
        return Advice(case, model_year, N_now, float(N_now_fit), float(N_target),
                      config, which, float(window), since, tuple(detected), ice,
                      temps or {}, factor, raw, clipped, N_after, skipped,
                      tuple(reasons), tuple(warnings), provenance)

    if not have_hi:
        reasons.append(f"no {ICE_VOLUME_VAR} series in the trend output")
        return result()
    if since is not None and n_usable < config.post_jump_min_years:
        reasons.append(f"only {n_usable} settled years since the jump at model year "
                       f"{since} (need {config.post_jump_min_years}, after "
                       f"{config.settle_years} settling years): keep running")
        return result()

    _, h_fit = _annual(columns, ICE_VOLUME_VAR, which)
    _, h_nat = _annual(columns, ICE_VOLUME_VAR, "native")
    h_now = float(h_nat[-1])

    # Fit once to locate the current point on the conduction law, then set the
    # target relative to it (the native N_now is interannually noisy).
    probe = extrapolate(ICE_VOLUME_VAR, N_fit[w], h_fit[w], N_now, h_now,
                        "hyperbolic", gate)
    if "hyperbolic" not in probe.fits:
        reasons.extend(probe.reasons)
        return result(ice=probe)
    a, b = probe.fits["hyperbolic"].params
    N_now_fit = a + b / h_now
    if config.N_target is not None:
        N_target = float(config.N_target)
        # An explicit target must advance toward equilibrium (the fitted
        # asymptote a) without crossing or reversing past the current point:
        # strictly between N_now_fit and a. The extrapolate() gate below
        # separately refuses a target at/past the asymptote (X -> inf); this
        # additionally catches a target on the correct side of the asymptote
        # but the wrong side of N_now_fit (moving further from equilibrium).
        lo_ok, hi_ok = sorted((N_now_fit, a))
        if not (lo_ok < N_target < hi_ok):
            reasons.append(f"explicit N_target {N_target:+.2f} does not lie strictly "
                           f"between the current imbalance N_now_fit={N_now_fit:+.2f} "
                           f"and the fitted asymptote a={a:+.2f}: refusing (would not "
                           f"advance toward equilibrium without crossing/reversing)")
    else:
        N_target = N_now_fit * (1.0 - config.n_fraction)
    # Clip the target to the trusted extrapolation range (the phase-space
    # analogue of the per-step magnitude clip; later jumps come out smaller,
    # like a decreasing dt schedule). Hindcast errors are flat up to 5x.
    Nw = N_fit[w]
    reach = gate.max_extrapolation_ratio * float(np.ptp(Nw))
    lo_edge, hi_edge = float(Nw.min()), float(Nw.max())
    if N_target > hi_edge + reach or N_target < lo_edge - reach:
        wanted = N_target
        N_target = hi_edge + reach if N_target > hi_edge else lo_edge - reach
        warnings.append(f"target N {wanted:+.2f} limited to {N_target:+.2f} by the "
                        f"window's N-range ({gate.max_extrapolation_ratio:g}x "
                        f"{np.ptp(Nw):.2f} W/m2)")
    ice = extrapolate(ICE_VOLUME_VAR, N_fit[w], h_fit[w], N_target, h_now,
                      "hyperbolic", gate)
    reasons.extend(ice.reasons)

    # Fail closed (feasibility-review finding 2): a missing ICEFRAC series
    # means the settled-edge assumption behind the conduction law and the
    # fixed-area vicen/eicen scaling cannot be confirmed at all, so refuse
    # rather than merely warn.
    if area is None:
        reasons.append(f"no {ICE_AREA_VAR} series: cannot confirm the ice edge "
                       f"settled — refusing (the conduction law and fixed-area "
                       f"scaling both assume a settled edge)")
    else:
        d_area = _icefrac_change(t, area, w)
        if d_area > config.max_icefrac_change:
            reasons.append(f"ice edge still moving: {ICE_AREA_VAR} changed by "
                           f"{d_area:.3f} across the window (> "
                           f"{config.max_icefrac_change}); the "
                           f"conduction law and fixed-area scaling need a settled edge")

    # Likewise a missing qi series, or qi/hi drifting more than the allowed
    # fraction, means enthalpy-per-volume stability (which scaling eicen with
    # vicen assumes) cannot be confirmed — refuse rather than warn.
    if not _have(columns, ICE_ENTHALPY_VAR, which):
        reasons.append(f"no {ICE_ENTHALPY_VAR} series: cannot confirm enthalpy per "
                       f"unit ice volume is steady before scaling eicen with vicen")
    else:
        _, q_fit = _annual(columns, ICE_ENTHALPY_VAR, which)
        ratio = q_fit[w] / h_fit[w]
        drift = float(abs(np.polyfit(t[w], ratio, 1)[0]) * window / abs(ratio.mean()))
        if drift > config.max_enthalpy_drift:
            reasons.append(f"{ICE_ENTHALPY_VAR}/{ICE_VOLUME_VAR} drifted {drift:.1%} "
                           f"over the window (> {config.max_enthalpy_drift:.0%}): "
                           f"enthalpy per unit ice volume is not steady, refusing to "
                           f"scale eicen with vicen")

    factor = raw = N_after = skipped = None
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
            skipped = float((h_new**2 - h_now**2) / h2_slope)

    temps = {}
    N_ref = N_after if N_after is not None else N_target
    for var in TEMPERATURE_VARS:
        if not _have(columns, var, which):
            continue
        _, X_fit = _annual(columns, var, which)
        _, X_nat = _annual(columns, var, "native")
        temps[var] = extrapolate(var, N_fit[w], X_fit[w], N_ref, float(X_nat[-1]),
                                 "linear", gate)
        if not temps[var].accepted:
            # A rejected temperature reference (kink / feedback threshold, or
            # too far beyond the sampled range) is fail-closed: refuse the
            # jump, do not merely warn (feasibility-review finding 2). The
            # `check` step relies on this reference being trustworthy.
            reasons.append(f"{var}(N) reference rejected: "
                           + "; ".join(temps[var].reasons))

    if any(not r.accepted for r in temps.values()):
        factor = raw = N_after = skipped = None
        clipped = False

    return result(N_target, N_now_fit, ice, temps, factor, raw, clipped, N_after,
                  skipped)
