"""Projection jump: carry a still-warming hot run to its trend-projected state
at a target model year in one coupled (somtp + cam.r) jump.

Pure computation over annual means (read via ``trend_io`` in the CLI). The
other advisors size a jump toward N = 0 and stay inside the hard bounds
(25 K, 50-yr probes); this one does not look for an equilibrium. It projects
the surface and the vapor column along their recent linear trends to
``target_year`` and builds the coupled state there (user request 2026-11-01,
atlasfu S5: "jumping 300 years is a gamechanger"). It therefore lifts the
bounds, and the jump must be confirmed with ``--override-limits``.

Why not each level's own trend. On S5 (years 126-145) the lower troposphere
warms ~2x faster than the surface and the top levels cool; extrapolated 255 yr
that gives a 503 K layer 64 K above a 435 K surface and 60 K at the top —
outside the k-tables, and a transient, not a state. So the column is
*shifted*:

* somtp and every level below ``taper_bottom_Pa`` warm by the projected dTS;
  the warming ramps to zero at ``ceiling_Pa`` (log p), the stratosphere and
  above are left alone — the same ceiling as the coupled jump (``atmos``).
* vapor: each layer's vapor mass per unit dry mass is multiplied by
  ``1 + (F - 1) w`` (F the projected / current vapor column, w the same
  ramp); the jump caps it at saturation at the new T, cell by cell.

Stability adjustments, each recorded in the advice:

* lapse rate — the taper cools the column's top relative to its bottom, which
  can make the ramp superadiabatic; the ceiling is raised (to lower pressure)
  by ``ceiling_step`` until the jumped mean column has no new superadiabatic
  layer (dry potential temperature, kappa = R/cp of the run);
* temperature range — the jumped mean column must stay ``T_margin`` inside
  the k-table range (``atmos.KTABLE_T_RANGE``); the jump re-checks every cell;
* vapor — saturation cap in the jump (above), previewed here on the mean
  column.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .atmos import KTABLE_T_RANGE, ceiling_weight, qsat

MODE = "project"


@dataclass(frozen=True)
class ProjectConfig:
    #: years of annual means the linear trends are fitted over, ending at the
    #: last year before the restart
    window: int = 20
    #: troposphere-only ceiling and full-strength bottom (Pa), as in ``atmos``
    ceiling_Pa: float = 1.0e4
    taper_bottom_Pa: float = 2.0e4
    #: margin kept inside the k-table range on the jumped mean column (K)
    T_margin: float = 10.0
    #: lapse-rate adjustment: divide the ceiling pressure by this until stable
    ceiling_step: float = 1.5
    min_ceiling_Pa: float = 100.0
    #: a level counts as superadiabatic when theta falls upward by more than this
    theta_tol: float = 0.1

    def __post_init__(self) -> None:
        if self.window < 5:
            raise ValueError("window must be >= 5 years")
        if not (0.0 < self.ceiling_Pa < self.taper_bottom_Pa):
            raise ValueError("need 0 < ceiling_Pa < taper_bottom_Pa")
        if self.ceiling_step <= 1.0:
            raise ValueError("ceiling_step must be > 1")


def linear_projection(years, values, end_year: int, target_year: int, window: int):
    """Least-squares line through ``values`` (n,) or (n, lev) over the
    ``window`` years ending at ``end_year``. Returns (slope, value at
    end_year, line at target_year); the 'now' value is the annual mean, not
    the fit."""
    years = np.asarray(years)
    v = np.asarray(values, dtype=float)
    if end_year not in years:
        raise ValueError(f"no annual mean for model year {end_year}")
    m = (years > end_year - window) & (years <= end_year)
    if m.sum() < window:
        raise ValueError(f"only {int(m.sum())} complete years in the {window}-yr "
                         f"window ending {end_year}")
    a = np.polyfit(years[m].astype(float), v[m], 1)
    now = v[list(years).index(end_year)]
    return a[0], now, a[0] * target_year + a[1]


def potential_temperature(T, p, kappa: float, p0: float = 1.0e5):
    return np.asarray(T, dtype=float) * (p0 / np.asarray(p, dtype=float)) ** kappa


def superadiabatic_levels(T, p, kappa: float, tol: float = 0.1) -> list:
    """Indices k (model order, top first) where theta(k) < theta(k+1) - tol:
    the layer from k+1 up to k is superadiabatic."""
    th = potential_temperature(T, p, kappa)
    return [int(k) for k in np.nonzero(th[:-1] < th[1:] - tol)[0]]


def project_state(years_TS, TS, years_TMQ, TMQ, years_lev, T_lev, Q_lev, P_lev,
                  end_year: int, target_year: int, kappa: float, epsilon: float,
                  config: ProjectConfig = ProjectConfig()) -> dict:
    """The shifted, vapor-scaled column at ``target_year`` (see module doc).

    ``T_lev``/``Q_lev``/``P_lev`` are (nyears, lev) global-mean annual means,
    model order (top first). Returns dT per level (``gain x dTS``), the vapor
    scale per level, and the bookkeeping for the advice."""
    c = config
    if target_year <= end_year:
        raise ValueError(f"target year {target_year} is not after {end_year}")
    sTS, TS_now, TS_tgt = linear_projection(years_TS, TS, end_year, target_year, c.window)
    sW, W_now, W_tgt = linear_projection(years_TMQ, TMQ, end_year, target_year, c.window)
    dTS = float(TS_tgt - TS_now)
    F = float(W_tgt / W_now)
    if F <= 0:
        raise ValueError(f"projected vapor column {W_tgt:.0f} kg/m2 is not positive")
    i = list(years_lev).index(end_year)
    T_now = np.asarray(T_lev, dtype=float)[i]
    Q_now = np.asarray(Q_lev, dtype=float)[i]
    p = np.asarray(P_lev, dtype=float)[i]

    before = set(superadiabatic_levels(T_now, p, kappa, c.theta_tol))
    ceiling = c.ceiling_Pa
    adjustments = []
    while True:
        w = ceiling_weight(p, ceiling, c.taper_bottom_Pa)
        T_new = T_now + dTS * w
        new = sorted(set(superadiabatic_levels(T_new, p, kappa, c.theta_tol)) - before)
        if not new:
            break
        nxt = ceiling / c.ceiling_step
        if nxt < c.min_ceiling_Pa:
            raise ValueError(f"no ceiling above {c.min_ceiling_Pa:g} Pa keeps the "
                             f"jumped column free of superadiabatic layers "
                             f"(levels {[k + 1 for k in new]})")
        adjustments.append(f"lapse rate: levels {[k + 1 for k in new]} superadiabatic "
                           f"with the ceiling at {ceiling / 100:.3g} hPa; raised to "
                           f"{nxt / 100:.3g} hPa")
        ceiling = nxt

    lo, hi = KTABLE_T_RANGE
    moved = T_new[w > 0]
    if moved.max() > hi - c.T_margin or moved.min() < lo + c.T_margin:
        raise ValueError(f"jumped levels of the mean column span {moved.min():.1f}-"
                         f"{moved.max():.1f} K, closer than {c.T_margin:g} K to the "
                         f"k-table range {lo:g}-{hi:g} K: choose an earlier target year")
    vapor_scale = 1.0 + (F - 1.0) * w
    r = Q_now / (1.0 - Q_now) * vapor_scale
    q_new = r / (1.0 + r)
    qs = qsat(T_new, p, epsilon)
    sat = np.nonzero(q_new > np.maximum(qs, Q_now))[0]
    if sat.size:
        adjustments.append(f"vapor: levels {[int(k) + 1 for k in sat]} would be "
                           f"supersaturated on the mean column; the jump caps q at "
                           f"saturation cell by cell")
    top = int(np.argmax(w >= 1.0)) if np.any(w >= 1.0) else int(p.size - 1)
    return {
        "dTS": dTS, "TS_now": float(TS_now), "TS_target": float(TS_tgt),
        "trend_K_per_yr": float(sTS),
        "TMQ_now": float(W_now), "TMQ_target": float(W_tgt),
        "TMQ_trend_per_yr": float(sW), "vapor_factor": F,
        "gain": w, "dT_levels": dTS * w, "vapor_scale": vapor_scale,
        "T_now": T_now, "T_after": T_new, "p_mid": p,
        "q_now": Q_now, "q_after": q_new,
        "ceiling_Pa": float(ceiling), "ceiling_Pa_requested": float(c.ceiling_Pa),
        "taper_bottom_Pa": float(c.taper_bottom_Pa), "top_index": top,
        "superadiabatic_before": sorted(int(k) + 1 for k in before),
        "adjustments": adjustments,
        "T_after_max": float(T_new.max()), "T_after_min": float(T_new.min()),
        "years": [int(end_year - c.window + 1), int(end_year)],
    }


def project_advice(case: str, end_year: int, target_year: int, res: dict,
                   config: ProjectConfig, schema_version: str,
                   provenance: Optional[dict] = None) -> dict:
    """Coupled advice (the ``atm-profile`` schema) for ``jump``; the
    ``atmosphere`` block is completed with the profile file by the CLI."""
    restart = f"{end_year + 1:04d}-01-01-00000"
    dTS = res["dTS"]
    return {
        "schema_version": schema_version, "plugin": "som_ocean", "mode": MODE,
        "case": case, "model_year": int(end_year), "since_year": None,
        "override_limits": True,
        "config": {"window": config.window, "target_year": int(target_year),
                   "heat_ratio": 1.0, "max_dT": None},
        "provenance": provenance or {},
        "TS_now": res["TS_now"], "trend_K_per_yr": res["trend_K_per_yr"],
        "TS_target": res["TS_target"], "TS_after": res["TS_now"] + dTS,
        "expected_after_jump": {"TS": res["TS_now"] + dTS},
        "TMQ_now": res["TMQ_now"], "TMQ_target": res["TMQ_target"],
        "vapor_factor": res["vapor_factor"],
        "somtp_dT_raw": dTS, "somtp_dT": dTS, "somtp_dT_ocean_only": None,
        "somtp_dT_clipped": False, "dTS": dTS, "sizing": "projection",
        "years_skipped": float(target_year - end_year),
        "reasons": [], "warnings": [
            f"projection jump: hard bounds lifted (somtp {dTS:+.1f} K); the target "
            f"is a linear trend extrapolated {target_year - end_year} yr, not an "
            f"equilibrium"] + list(res["adjustments"]),
        "atmosphere": {
            "method": MODE, "years": res["years"], "restart_date": restart,
            "dTS": dTS, "levels": int(res["gain"].size),
            "top_level": res["top_index"],
            "top_pressure_Pa": float(res["p_mid"][res["top_index"]]),
            "ceiling_Pa": res["ceiling_Pa"],
            "ceiling_Pa_requested": res["ceiling_Pa_requested"],
            "taper_bottom_Pa": res["taper_bottom_Pa"],
            "humidity": "vapor scale (capped at saturation)",
            "vapor_scale_surface": float(res["vapor_scale"][-1]),
            "T_after_max": res["T_after_max"], "T_after_min": res["T_after_min"],
            "superadiabatic_before": res["superadiabatic_before"],
            "adjustments": list(res["adjustments"]),
        },
    }
