"""cam_atm plugin: shift the atmosphere's temperature and water vapor in cam.r.

The companion of ``som_ocean`` for hot, thick, moist atmospheres. Pure arrays;
the netCDF read/write lives in ``restart``.

Why the atmosphere must move with the ocean (docs/ocean-jump.md, Finding 1):
at 365-375 K the atmosphere's heat capacity — almost all of it the latent
heat of the water-vapor column — is 2-3x the 50 m slab's. A somtp-only jump
is drained into the atmosphere within weeks and lands at C_ocean/C_total of
its size. Warming T alone is not enough either: at fixed q the air is left
undersaturated and the ocean re-evaporates to refill it. So T and q move
together, q at fixed relative humidity.

What cam.r holds, and what the edit must keep consistent (CAM FV, ExoCAM
cesm1.2.1; verified on atlasfu D4 0051):

* ``PT`` = T_v / pkz — scaled virtual potential temperature, with
  ``T_v = T (1 + zvir q)`` and ``pkz`` the layer mean of p^kappa
  (``dyn_comp.F90``, ``te_map.F90``); pressures from ``DELP`` down from ``ptop``;
* ``Q`` is moist specific humidity; ``DELP`` is moist. CAM's own
  ``dme_adjust`` adds mass when vapor is added, so the jump keeps each layer's
  *dry* mass: ``DELP' = DELP (1 - q) / (1 - q')``; ``PS`` follows. Condensate
  mixing ratios (``CLDLIQ``, ``CLDICE``, ``LCWAT``) are diluted so the
  condensate mass is unchanged;
* ``TEOUT`` — the column total energy at the end of the last physics step.
  The global energy fixer (``check_energy_gmean``) heats or cools the whole
  atmosphere by (energy now - TEOUT) on the first step, so a jump that left it
  alone would be undone at once. It is recomputed with CAM's own formula
  (``total_energy``: dry static energy cp T + g z + phis, kinetic energy,
  (latvap+latice) q, latice liq), which reproduces the stored TEOUT to 0.005 %;
* ``TCWAT`` / ``QCWAT`` (previous-step T and q of the stratiform scheme) and
  ``T_TTEND`` shift with T and q so the first step sees no spurious forcing;
* ``U``, ``V`` and ``PHIS`` are unchanged.

The increment's vertical shape is measured (``measured_profile``): the
horizontal-mean warming per level between two archived cam.i, per K of
surface warming (the per-layer horizontal-mean design this tool started
from). It is applied from the surface up to where the measured warming first
turns negative (the stratosphere cools; it holds negligible heat and adjusts
radiatively), with a linear taper over ``taper_levels`` above that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np

from .plugins import PLUGIN_REGISTRY, ConstraintReport, VariablePlugin

#: CAM4 physconst latent heats (shr_const), J/kg
LATVAP = 2.501e6
LATICE = 3.337e5
#: water-vapor gas constant, J/kg/K
RH2O = 461.504639820160

#: Fields the jump reads from and writes to cam.r.
STATE_FIELDS = ("PT", "Q", "DELP", "PS", "U", "V", "PHIS", "CLDLIQ", "CLDICE",
                "TEOUT", "TCWAT", "QCWAT", "LCWAT", "T_TTEND")
WRITTEN_FIELDS = ("PT", "Q", "DELP", "PS", "CLDLIQ", "CLDICE", "TEOUT", "TCWAT",
                  "QCWAT", "LCWAT", "T_TTEND")

#: Hard bound on any level's temperature increment, K.
MAX_LEVEL_DT = 25.0
#: Where saturation vapor pressure exceeds this fraction of the pressure the
#: fixed-RH formula is not trusted (near boiling): q is left alone there.
MAX_ES_FRACTION = 0.9
#: Refuse when the constants reproduce the stored TEOUT worse than this.
MAX_TE_MISMATCH = 1e-3
#: Refuse a jump that raises the global-mean vapor column by more than this
#: fraction (a +4.9 K coupled probe on atlasfu D4 adds 34 %).
MAX_VAPOR_INCREASE = 0.5


@dataclass(frozen=True)
class AtmConstants:
    """Model constants (from the run's atm.log: CPDAIR, RAIR, ZVIR, PTOP,
    SURFACE GRAVITY)."""
    cpair: float
    rair: float
    zvir: float
    gravit: float
    ptop: float
    latvap: float = LATVAP
    latice: float = LATICE

    @property
    def kappa(self) -> float:
        return self.rair / self.cpair

    @property
    def epsilon(self) -> float:
        """Rd / Rv — the q_s(e_s, p) molecular-weight ratio."""
        return 1.0 / (1.0 + self.zvir)


# ---------------------------------------------------------------------------
# column thermodynamics (arrays shaped (lev, ...) with lev top-down)
# ---------------------------------------------------------------------------

def interfaces(delp, ptop: float) -> np.ndarray:
    delp = np.asarray(delp, dtype=float)
    return ptop + np.concatenate([np.zeros((1,) + delp.shape[1:]), np.cumsum(delp, 0)])


def pkz(pe, kappa: float) -> np.ndarray:
    """FV layer-mean p^kappa: (pe_k+1^k - pe_k^k) / (k ln(pe_k+1/pe_k))."""
    pk = pe ** kappa
    lnp = np.log(pe)
    return (pk[1:] - pk[:-1]) / (kappa * (lnp[1:] - lnp[:-1]))


def p_mid(pe) -> np.ndarray:
    return 0.5 * (pe[1:] + pe[:-1])


def temperature(PT, Q, delp, c: AtmConstants) -> np.ndarray:
    return np.asarray(PT) * pkz(interfaces(delp, c.ptop), c.kappa) / (1.0 + c.zvir * np.asarray(Q))


def potential_temperature_field(T, Q, delp, c: AtmConstants) -> np.ndarray:
    """PT for a given T, q and layer thickness."""
    return np.asarray(T) * (1.0 + c.zvir * np.asarray(Q)) / pkz(interfaces(delp, c.ptop), c.kappa)


def mid_heights(T, Q, delp, c: AtmConstants) -> np.ndarray:
    """Layer-midpoint heights above the surface (CAM geopotential_t)."""
    pe = interfaces(delp, c.ptop)
    Tv = np.asarray(T) * (1.0 + c.zvir * np.asarray(Q))
    lnp = np.log(pe)
    hkl = lnp[1:] - lnp[:-1]
    hkk = 1.0 - pe[:-1] * hkl / delp
    rog = c.rair / c.gravit
    zm = np.empty_like(Tv)
    zi = np.zeros_like(pe[0])
    for k in range(Tv.shape[0] - 1, -1, -1):
        zm[k] = zi + rog * Tv[k] * hkk[k]
        zi = zi + rog * Tv[k] * hkl[k]
    return zm


def total_energy(T, Q, delp, U, V, PHIS, CLDLIQ, c: AtmConstants) -> np.ndarray:
    """Column total energy exactly as CAM's check_energy (J/m2), per column."""
    zm = mid_heights(T, Q, delp, c)
    s = c.cpair * T + c.gravit * zm + PHIS
    g = c.gravit
    return (((s + 0.5 * (U ** 2 + V ** 2)) * delp).sum(0) / g
            + (c.latvap + c.latice) * (Q * delp).sum(0) / g
            + c.latice * (CLDLIQ * delp).sum(0) / g)


def esat(T) -> np.ndarray:
    """Saturation vapor pressure (Pa): Wagner & Pruss (2002) over liquid above
    273.16 K, Murphy & Koop (2005) over ice below. Only ratios e_s(T')/e_s(T)
    enter the fixed-RH update, so the absolute calibration matters little."""
    T = np.asarray(T, dtype=float)
    Tc, Pc = 647.096, 22.064e6
    tau = 1.0 - np.minimum(T, Tc - 1e-6) / Tc
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    ln_liq = (Tc / np.minimum(T, Tc - 1e-6)) * (
        a[0] * tau + a[1] * tau ** 1.5 + a[2] * tau ** 3 + a[3] * tau ** 3.5
        + a[4] * tau ** 4 + a[5] * tau ** 7.5)
    liq = Pc * np.exp(ln_liq)
    Ti = np.maximum(T, 1.0)
    ice = np.exp(9.550426 - 5723.265 / Ti + 3.53068 * np.log(Ti) - 0.00728332 * Ti)
    return np.where(T >= 273.16, liq, ice)


def qsat(T, p, epsilon: float) -> np.ndarray:
    es = esat(T)
    return epsilon * es / np.maximum(p - (1.0 - epsilon) * es, 1e-3 * p)


# ---------------------------------------------------------------------------
# the measured vertical profile
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ProfileConfig:
    #: 3-point vertical smoothing passes
    smooth: int = 1
    #: bound on the per-level warming per K of surface warming
    max_gain: float = 4.0
    #: "auto" = up to where the measured gain first turns negative from the
    #: surface; or a pressure (Pa) above which nothing is changed
    top: object = "auto"
    #: linear taper to zero over this many levels above the top
    taper_levels: int = 2
    #: troposphere-only ceiling (Pa): nothing is changed at pressures below
    #: ``ceiling_Pa`` (the stratosphere and above are numerically fragile), and
    #: the jump ramps from full strength at ``taper_bottom_Pa`` to zero at the
    #: ceiling, linearly in log p. Applies on top of ``top``.
    ceiling_Pa: float = 1.0e4
    taper_bottom_Pa: float = 2.0e4
    #: refuse a profile measured over less surface warming than this (K)
    min_dTS: float = 0.5

    def __post_init__(self) -> None:
        if self.smooth < 0 or self.taper_levels < 0:
            raise ValueError("smooth and taper_levels must be >= 0")
        if not (0.0 < self.max_gain <= 10.0):
            raise ValueError("max_gain must be in (0, 10]")
        if self.top != "auto" and not (isinstance(self.top, (int, float)) and self.top > 0):
            raise ValueError(f"top must be 'auto' or a pressure in Pa, got {self.top!r}")
        if not (0.0 < self.ceiling_Pa < self.taper_bottom_Pa):
            raise ValueError("need 0 < ceiling_Pa < taper_bottom_Pa")


def ceiling_weight(p, ceiling_Pa: float, taper_bottom_Pa: float) -> np.ndarray:
    """1 at p >= taper_bottom_Pa, 0 at p <= ceiling_Pa, linear in log p between."""
    p = np.asarray(p, dtype=float)
    x = (np.log(p) - np.log(ceiling_Pa)) / (np.log(taper_bottom_Pa) - np.log(ceiling_Pa))
    return np.clip(x, 0.0, 1.0)


@dataclass(frozen=True)
class AtmProfile:
    gain: np.ndarray          # (lev,) K of warming per K of surface warming, applied
    raw_gain: np.ndarray      # (lev,) measured, before smoothing/clipping/taper
    p_mid: np.ndarray         # (lev,) global-mean layer pressure, Pa
    T_now: np.ndarray         # (lev,) horizontal-mean T at the later restart
    dTS: float                # surface warming the gain is normalized by
    dt_years: float
    top_index: int            # highest level (smallest index) with full weight
    config: ProfileConfig

    def summary(self) -> dict:
        k0 = self.top_index
        return {"dTS_baseline": self.dTS, "baseline_years": self.dt_years,
                "levels": int(self.gain.size), "top_level": int(k0),
                "top_pressure_Pa": float(self.p_mid[k0]),
                "gain_surface": float(self.gain[-1]),
                "gain_max": float(self.gain.max()),
                "gain_mass_weighted": float(np.average(self.gain, weights=np.gradient(self.p_mid))),
                "smooth": self.config.smooth, "max_gain": self.config.max_gain,
                "top": self.config.top, "taper_levels": self.config.taper_levels,
                "ceiling_Pa": self.config.ceiling_Pa,
                "taper_bottom_Pa": self.config.taper_bottom_Pa,
                "gain_above_ceiling_max": float(
                    np.abs(self.gain[self.p_mid < self.config.ceiling_Pa]).max(initial=0.0))}


def measured_profile(T_old, T_new, dTS: float, p_mid_now, dt_years: float,
                     config: ProfileConfig = ProfileConfig()) -> AtmProfile:
    """Gain per level from two horizontal-mean temperature profiles.

    ``T_old``/``T_new`` are (lev,) area-weighted means at the two restarts,
    ``dTS`` the area-mean surface (somtp) warming between them.
    """
    T_old = np.asarray(T_old, dtype=float)
    T_new = np.asarray(T_new, dtype=float)
    p = np.asarray(p_mid_now, dtype=float)
    if T_old.shape != T_new.shape or T_new.ndim != 1 or p.shape != T_new.shape:
        raise ValueError("profiles must be matching 1-D (lev,) arrays")
    if not (np.all(np.isfinite(T_old)) and np.all(np.isfinite(T_new))):
        raise ValueError("non-finite temperatures in the profiles")
    if abs(dTS) < config.min_dTS:
        raise ValueError(f"surface warmed only {dTS:+.2f} K between the restarts "
                         f"(< {config.min_dTS:g}): the profile would be noise — use "
                         f"restarts further apart")
    raw = (T_new - T_old) / dTS
    g = raw.copy()
    for _ in range(config.smooth):
        g = np.r_[g[0], (g[:-2] + 2 * g[1:-1] + g[2:]) / 4.0, g[-1]]
    n = g.size
    if config.top == "auto":
        k0 = n - 1
        while k0 > 0 and g[k0 - 1] > 0:
            k0 -= 1
    else:
        above = np.where(p >= float(config.top))[0]
        k0 = int(above[0]) if above.size else n - 1
    w = np.zeros(n)
    w[k0:] = 1.0
    for j in range(1, config.taper_levels + 1):
        if k0 - j >= 0:
            w[k0 - j] = 1.0 - j / (config.taper_levels + 1)
    g = np.clip(g, 0.0, config.max_gain)
    # above the top the measured gain is negative (stratosphere) or noise: ramp
    # the top level's gain down to zero instead
    g = np.where(np.arange(n) < k0, g[k0] * w, g)
    # troposphere only: zero at and above the ceiling, log-p ramp below it
    wc = ceiling_weight(p, config.ceiling_Pa, config.taper_bottom_Pa)
    g = g * wc
    full = np.where((wc >= 1.0) & (np.arange(n) >= k0))[0]
    k_full = int(full[0]) if full.size else n - 1
    return AtmProfile(g, raw, p, T_new, float(dTS), float(dt_years), k_full, config)


# ---------------------------------------------------------------------------
# the plugin
# ---------------------------------------------------------------------------

def jump_state(state: Dict[str, np.ndarray], dT_levels, c: AtmConstants,
               iterations: int = 4) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """Raise T by ``dT_levels`` (per level, broadcast over columns; K) and q at
    fixed relative humidity, keeping each layer's dry mass. Returns the new
    written fields and a report."""
    PT, Q, D = (np.asarray(state[k], dtype=float) for k in ("PT", "Q", "DELP"))
    dT = np.asarray(dT_levels, dtype=float)
    if dT.ndim == 1:
        if dT.size != PT.shape[0]:
            raise ValueError(f"increment has {dT.size} levels, cam.r has {PT.shape[0]}")
        dT = dT.reshape((-1,) + (1,) * (PT.ndim - 1))
    if not np.all(np.isfinite(dT)) or np.max(np.abs(dT)) > MAX_LEVEL_DT:
        raise ValueError(f"level increment outside +/-{MAX_LEVEL_DT:g} K or non-finite")
    for k, v in state.items():
        if not np.all(np.isfinite(v)):
            raise ValueError(f"cam.r {k} has non-finite values")
    if np.any(Q < 0) or np.any(Q >= 1) or np.any(D <= 0):
        raise ValueError("cam.r Q or DELP out of range")

    T = temperature(PT, Q, D, c)
    Tn = T + dT
    pe = interfaces(D, c.ptop)
    pm = p_mid(pe)
    es0 = esat(T)
    near_boil = (es0 > MAX_ES_FRACTION * pm) | (esat(Tn) > MAX_ES_FRACTION * pm)
    dry = D * (1.0 - Q)
    qs0 = qsat(T, pm, c.epsilon)
    Qn, Dn = Q.copy(), D.copy()
    for _ in range(iterations):
        pmn = p_mid(interfaces(Dn, c.ptop))
        Qn = np.where(near_boil, Q, Q * qsat(Tn, pmn, c.epsilon) / qs0)
        Qn = np.minimum(Qn, 0.95)
        Dn = dry / (1.0 - Qn)
    PTn = potential_temperature_field(Tn, Qn, Dn, c)
    dil = D / Dn
    te0 = total_energy(T, Q, D, state["U"], state["V"], state["PHIS"], state["CLDLIQ"], c)
    ten = total_energy(Tn, Qn, Dn, state["U"], state["V"], state["PHIS"],
                       state["CLDLIQ"] * dil, c)
    teout = np.asarray(state["TEOUT"], dtype=float)
    mismatch = float(np.abs(te0 / teout - 1.0).max())
    out = {
        "PT": PTn, "Q": Qn, "DELP": Dn,
        "PS": interfaces(Dn, c.ptop)[-1],
        "CLDLIQ": np.asarray(state["CLDLIQ"]) * dil,
        "CLDICE": np.asarray(state["CLDICE"]) * dil,
        "LCWAT": np.asarray(state["LCWAT"]) * dil,
        "TEOUT": teout + (ten - te0),
        "TCWAT": np.asarray(state["TCWAT"]) + (Tn - T),
        "T_TTEND": np.asarray(state["T_TTEND"]) + (Tn - T),
        "QCWAT": np.asarray(state["QCWAT"]) + (Qn - Q),
    }
    g = c.gravit
    report = {
        "te_mismatch_max": mismatch,
        "dTE_mean_J_m2": float((ten - te0).mean()),
        "vapor_added_kg_m2": float(((Qn * Dn - Q * D).sum(0) / g).mean()),
        "vapor_increase_fraction": float(((Qn * Dn - Q * D).sum(0)).mean()
                                         / ((Q * D).sum(0)).mean()),
        "dPS_mean_Pa": float((out["PS"] - np.asarray(state["PS"])).mean()),
        "dT_max": float(np.max(Tn - T)), "dT_min": float(np.min(Tn - T)),
        "near_boiling_cells": int(near_boil.sum()),
        "q_capped_cells": int((Qn >= 0.95).sum()),
        "T_check_rms": float(np.sqrt(np.mean(
            (temperature(PTn, Qn, Dn, c) - Tn) ** 2))),
    }
    return out, report


class CamAtmPlugin(VariablePlugin):
    name = "cam_atm"
    file_kind = "cam.r"
    variables = WRITTEN_FIELDS

    def __init__(self, constants: Optional[AtmConstants] = None):
        self.constants = constants

    def apply_delta(self, fields, delta):
        if self.constants is None:
            raise ValueError("cam_atm needs the model constants (AtmConstants)")
        out, self.last_report = jump_state(fields, delta, self.constants)
        return out

    def enforce_constraints(self, fields_before, fields_after):
        """The physics is enforced inside ``jump_state`` (dry mass per layer,
        fixed RH, energy bookkeeping); this pass refuses a result whose
        constants did not reproduce the model's own stored energy."""
        rep = getattr(self, "last_report", {})
        if rep.get("te_mismatch_max", 0.0) > MAX_TE_MISMATCH:
            raise ValueError(f"the model constants reproduce the stored TEOUT only to "
                             f"{rep['te_mismatch_max']:.2%} (> {MAX_TE_MISMATCH:.1%}): "
                             f"wrong constants for this case — refusing")
        if rep.get("vapor_increase_fraction", 0.0) > MAX_VAPOR_INCREASE:
            raise ValueError(f"the jump would add {rep['vapor_increase_fraction']:.0%} "
                             f"to the vapor column (> {MAX_VAPOR_INCREASE:.0%}): too "
                             f"large a shock — use a smaller surface step")
        adj = {k: (f"mean {float(np.mean(fields_before[k])):.6g} -> "
                   f"{float(np.mean(v)):.6g}") for k, v in fields_after.items()}
        return fields_after, ConstraintReport(plugin=self.name, adjustments=adj)


PLUGIN_REGISTRY[CamAtmPlugin.name] = CamAtmPlugin()
