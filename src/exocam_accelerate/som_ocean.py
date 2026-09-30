"""som_ocean plugin: shift the slab-ocean temperature ``somtp`` in docn.r.

The hot / thick-atmosphere counterpart of ``aqua_ice``. Pure arrays; the
netCDF read/write lives in ``restart``.

What docn.r holds (ExoCAM ``docn_comp_mod.F90``, SOM mode):

* one variable, ``somtp(gsize)`` — the prognostic mixed-layer temperature in
  **Kelvin** (the SSTDATA/SOM paths add ``TkFrz``);
* ``gsize = ni * nj`` of the docn domain, in the MCT global-index order
  ``n = i + (j - 1) * ni`` (lon fastest), so ``somtp.reshape(nj, ni)`` is the
  lat-lon map (``to_grid``);
* on a restart the first coupling step copies ``somtp`` straight into the
  ocean->coupler ``So_t`` (no update on the first call), so an edited somtp is
  the surface temperature the atmosphere sees from the first step on;
* land/masked cells (``imask == 0``) are never updated by the model; changing
  them is harmless but pointless, so they are left alone when a mask is known.

What is changed, and why this is self-consistent for a continuation:

* only ``somtp``: ``pop_frc`` (``hblt``, ``qdp``) is a forcing stream, and the
  stream restart ``docn.rs1.bin`` holds only stream-time bookkeeping;
* cells at the sea-water freezing point may carry sea ice in cice.r — warming
  them would put open-ocean heat under ice, so they are never changed (the
  advisor separately refuses cases with ice at all);
* the coupler snapshot (cpl.r ``o2x_So_t``, cam.rs ``x2a_Sx_t``) is stale for
  one coupling interval, like the ice fraction in the aqua_ice argument.

The increment ``dT`` is one number (uniform) or a per-cell map (the measured
warming pattern, ``pattern_weights``); either way its area mean over open
ocean is the ``somtp_dT`` the advisor sized.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from .plugins import PLUGIN_REGISTRY, ConstraintReport, VariablePlugin

SOM_FIELDS = ("somtp",)

#: shr_const_TkFrzSw: freezing point of sea water, K (TkFrz - 1.8)
TK_FRZ_SW = 273.15 - 1.8
#: cells within this of the freezing point are treated as (possibly) ice-covered
FREEZE_MARGIN = 0.05
#: Hard bound on any cell's increment, whatever the advisor says (K).
MAX_DT = 25.0
#: Plausible somtp range, K — outside it the file is not a SOM temperature.
SOMTP_RANGE = (150.0, 500.0)


def check_dT(dT, what: str = "somtp increment") -> np.ndarray:
    """A finite increment (scalar or map) within the hard bound."""
    d = np.asarray(dT, dtype=float)
    if not np.all(np.isfinite(d)):
        raise ValueError(f"{what} has non-finite values")
    if np.max(np.abs(d)) > MAX_DT:
        raise ValueError(f"{what} {np.max(np.abs(d)):g} K outside the hard bound "
                         f"+/-{MAX_DT:g} K")
    return d


def to_grid(vec, nj: int, ni: int) -> np.ndarray:
    """docn gsize vector -> (nj, ni) lat-lon map (lon fastest)."""
    v = np.asarray(vec)
    if v.shape != (nj * ni,):
        raise ValueError(f"somtp has shape {v.shape}, the domain is {nj}x{ni} "
                         f"= {nj * ni} cells")
    return v.reshape(nj, ni)


def area_mean(x, area, sel=None) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(area, dtype=float)
    if sel is not None:
        w = np.where(sel, w, 0.0)
    if w.sum() <= 0:
        raise ValueError("no cells with positive area to average over")
    return float((x * w).sum() / w.sum())


class SomOceanPlugin(VariablePlugin):
    name = "som_ocean"
    file_kind = "docn.r"
    variables = SOM_FIELDS

    def apply_delta(self, fields: Dict[str, np.ndarray], delta) -> Dict[str, np.ndarray]:
        """``delta`` is the increment in K: a scalar or an array shaped like somtp."""
        somtp = np.asarray(fields["somtp"], dtype=float)
        d = check_dT(delta)
        if d.ndim and d.shape != somtp.shape:
            raise ValueError(f"increment map is {d.shape}, somtp is {somtp.shape}")
        return {"somtp": somtp + d}

    def enforce_constraints(self, fields_before, fields_after, ocean=None):
        """Leave near-freezing cells and non-ocean cells unchanged; keep the rest
        above freezing. ``ocean`` (bool, somtp-shaped) marks cells the model
        updates; None = all.
        """
        before = np.asarray(fields_before["somtp"], dtype=float)
        after = np.array(fields_after["somtp"], dtype=float)
        if not np.all(np.isfinite(before)):
            raise ValueError("somtp has non-finite values in the restart")
        lo, hi = SOMTP_RANGE
        if before.min() < lo or before.max() > hi:
            raise ValueError(f"somtp spans {before.min():.1f}-{before.max():.1f}, "
                             f"outside {lo:g}-{hi:g} K: not a SOM temperature in K")
        frozen = before <= TK_FRZ_SW + FREEZE_MARGIN
        keep = frozen.copy()
        if ocean is not None:
            keep |= ~np.asarray(ocean, dtype=bool)
        after[keep] = before[keep]
        below = after < TK_FRZ_SW
        after[below] = TK_FRZ_SW
        d = after - before
        msg = (f"mean {before.mean():.3f} -> {after.mean():.3f} K "
               f"(increment {d.min():+.3f}..{d.max():+.3f} K)")
        if frozen.any():
            msg += f"; left {int(frozen.sum())} near-freezing cells unchanged"
        if ocean is not None and (~np.asarray(ocean, dtype=bool)).any():
            msg += f"; left {int((~np.asarray(ocean, dtype=bool)).sum())} non-ocean cells unchanged"
        if below.any():
            msg += f"; clamped {int(below.sum())} cells at the freezing point"
        return {"somtp": after}, ConstraintReport(plugin=self.name,
                                                  adjustments={"somtp": msg})


PLUGIN_REGISTRY[SomOceanPlugin.name] = SomOceanPlugin()


# ---------------------------------------------------------------------------
# measured warming pattern (optional; the default jump is uniform)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PatternConfig:
    #: box-smoothing half-width in cells (lon periodic); 0 = none. Jan-1
    #: snapshots carry weather noise that a 3x3 box average damps.
    smooth: int = 1
    #: per-cell weight bounds (area mean is 1)
    min_weight: float = 0.0
    max_weight: float = 3.0
    #: refuse when the area-mean rate is smaller than this (K/yr): the
    #: pattern would be mostly noise
    min_mean_rate: float = 0.02

    def __post_init__(self) -> None:
        if self.smooth < 0:
            raise ValueError(f"smooth must be >= 0, got {self.smooth!r}")
        if not (0.0 <= self.min_weight < 1.0 < self.max_weight):
            raise ValueError(f"need 0 <= min_weight < 1 < max_weight, got "
                             f"({self.min_weight!r}, {self.max_weight!r})")
        if self.min_mean_rate <= 0:
            raise ValueError("min_mean_rate must be positive")


@dataclass(frozen=True)
class OceanPattern:
    weight: np.ndarray       # (nj, ni), area mean 1 over open ocean; 0 elsewhere
    rate: np.ndarray         # smoothed d(somtp)/dt, K/yr
    somtp_now: np.ndarray    # K
    area: np.ndarray         # 0 off the open-ocean mask
    mean_rate: float
    dt_years: float
    config: PatternConfig

    def summary(self) -> dict:
        sel = self.area > 0
        w = self.weight[sel]
        a = self.area[sel] / self.area[sel].sum()
        return {"baseline_years": self.dt_years,
                "mean_rate_K_per_yr": self.mean_rate,
                "smooth": self.config.smooth,
                "weight_bounds": [self.config.min_weight, self.config.max_weight],
                "weight_min": float(w.min()), "weight_max": float(w.max()),
                "weight_std": float(np.sqrt((a * (w - 1.0) ** 2).sum())),
                "area_at_min": float(a[w <= self.config.min_weight + 1e-9].sum()),
                "area_at_max": float(a[w >= self.config.max_weight - 1e-9].sum())}


def _box_smooth(x, valid, k: int):
    """Mean over a (2k+1)^2 box of valid cells; lon periodic, lat clipped."""
    if k == 0:
        return np.where(valid, x, 0.0)
    nj, ni = x.shape
    v = valid.astype(float)
    xv = np.where(valid, x, 0.0)
    num = np.zeros_like(xv)
    den = np.zeros_like(xv)
    for dj in range(-k, k + 1):
        rows = np.clip(np.arange(nj) + dj, 0, nj - 1)
        for di in range(-k, k + 1):
            num += np.roll(xv[rows], -di, axis=1)
            den += np.roll(v[rows], -di, axis=1)
    return np.where(valid & (den > 0), num / np.maximum(den, 1.0), 0.0)


def pattern_weights(somtp_old, somtp_now, dt_years: float, area, ocean=None,
                    config: PatternConfig = PatternConfig()) -> OceanPattern:
    """Per-cell weights from the local warming rate between two restarts.

    ``weight = rate / area_mean(rate)``, clipped to the configured bounds and
    renormalized to area mean 1 over open ocean, so ``dT_cell = dT * weight``
    keeps the area-mean increment the advisor sized. Near-freezing and
    non-ocean cells get weight 0. Refuses (ValueError) when the mean rate is
    too small for the pattern to be more than noise.
    """
    old = np.asarray(somtp_old, dtype=float)
    now = np.asarray(somtp_now, dtype=float)
    area = np.asarray(area, dtype=float)
    if old.shape != now.shape or now.shape != area.shape or now.ndim != 2:
        raise ValueError(f"shapes differ or are not 2-D: {old.shape}, {now.shape}, "
                         f"{area.shape}")
    if dt_years <= 0:
        raise ValueError(f"dt_years must be positive, got {dt_years!r}")
    if not (np.all(np.isfinite(old)) and np.all(np.isfinite(now))):
        raise ValueError("somtp has non-finite values")
    sel = (area > 0) & (now > TK_FRZ_SW + FREEZE_MARGIN) & (old > TK_FRZ_SW + FREEZE_MARGIN)
    if ocean is not None:
        sel &= np.asarray(ocean, dtype=bool)
    if not sel.any():
        raise ValueError("no open-ocean cells to build a pattern on")
    rate = _box_smooth((now - old) / dt_years, sel, config.smooth)
    mean = area_mean(rate, area, sel)
    if abs(mean) < config.min_mean_rate:
        raise ValueError(f"area-mean somtp rate {mean:+.3f} K/yr is below "
                         f"{config.min_mean_rate:g}: the pattern would be noise — "
                         f"use the uniform jump")
    w = np.where(sel, rate / mean, 0.0)
    for _ in range(50):                  # clip, renormalize, repeat to a fixed point
        w = np.where(sel, np.clip(w, config.min_weight, config.max_weight), 0.0)
        m = area_mean(w, area, sel)
        if abs(m - 1.0) < 1e-10:
            break
        w = w / m
    area_used = np.where(sel, area, 0.0)
    return OceanPattern(w, rate, now, area_used, mean, float(dt_years), config)
