"""Tapered ice jump: scale only the ice that is still conduction-limited.

Pure arrays; the netCDF reads live in ``restart``. Decided 2026-09-29 from the
first two real jumps (grp4 pt01/pt03, x1.5 at 0141): the uniform factor is
right for the thick, Stefan-growing ice but overshoots ice that was already
near its local equilibrium (the substellar thin-ice ring) — there the extra
ice melts back instead of buying time.

Test for "conduction-limited", per grid cell, from two pristine restarts
``dt`` years apart (grid-box ice volume per area h = sum over categories of
``vicen``):

    S = dh/dt * h_mid                   Stefan: dh/dt = k dT / (rho L h), so S
                                        is ~constant for conduction-limited ice
    w = clip(S / S_ref, 0, 1)           S_ref = median S over the thicker half
                                        of fully ice-covered cells
    weight = ramp(w; w_lo, w_hi)        0 below w_lo, 1 above w_hi, linear between
    factor_cell = 1 + (F - 1) * weight  F = the peak factor

Cells with partial ice cover get weight 0 (fixed-area scaling needs full
cover). Retro-test on pt01/pt03 (data before the jump only): cells with
w < 0.05 melted back or stalled after the uniform x1.5; growth rose smoothly
to the Stefan rate by w ~ 0.4 — hence the default ramp (0.05, 0.40).

Sizing F for a target imbalance: total basal freezing carries the TOA deficit
(energy_bot ~ energy_top), and a cell's freezing rate scales as 1/h, so after
the jump

    (N_after - a) / (N_now - a) = R(F) = sum(A g / f) / sum(A g)

with g the pre-jump growth rate per cell, f its factor, and a the fitted
conduction-law asymptote. For a uniform factor R = 1/F, which is the global
law a + b/(F h) the advisor uses — the taper generalizes it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

#: Minimum ice area fraction for a cell to be scaled at all.
FULL_COVER = 0.99


@dataclass(frozen=True)
class TaperConfig:
    ramp_lo: float = 0.05
    ramp_hi: float = 0.40
    full_cover: float = FULL_COVER

    def __post_init__(self) -> None:
        if not (0.0 <= self.ramp_lo < self.ramp_hi <= 1.0):
            raise ValueError(f"need 0 <= ramp_lo < ramp_hi <= 1, got "
                             f"({self.ramp_lo!r}, {self.ramp_hi!r})")
        if not (0.0 < self.full_cover <= 1.0):
            raise ValueError(f"full_cover must be in (0, 1], got {self.full_cover!r}")


@dataclass(frozen=True)
class TaperMask:
    weight: np.ndarray        # (nj, ni) in [0, 1]; 0 = never scaled
    w: np.ndarray             # raw S / S_ref, clipped to [0, 1]
    growth: np.ndarray        # pre-jump dh/dt, m/yr
    h_now: np.ndarray         # grid-box ice volume per area at the jump, m
    area: np.ndarray          # cell area (0 off the ice/ocean mask)
    S_ref: float              # m^2/yr
    dt_years: float
    config: TaperConfig

    def area_fraction(self, sel) -> float:
        return float(self.area[sel].sum() / self.area.sum())

    def summary(self) -> dict:
        return {
            "S_ref_m2_per_yr": self.S_ref,
            "baseline_years": self.dt_years,
            "ramp": [self.config.ramp_lo, self.config.ramp_hi],
            "full_cover": self.config.full_cover,
            "area_full_weight": self.area_fraction(self.weight >= 1.0),
            "area_partial_weight": self.area_fraction((self.weight > 0) & (self.weight < 1)),
            "area_unscaled": self.area_fraction((self.weight <= 0) & (self.area > 0)),
            "mean_weight": float((self.weight * self.area).sum() / self.area.sum()),
        }


def stefan_mask(h_old, h_now, dt_years: float, aice, area,
                config: TaperConfig = TaperConfig()) -> TaperMask:
    """Per-cell taper weights from two pristine states ``dt_years`` apart."""
    h_old, h_now, aice, area = (np.asarray(x, dtype=float)
                                for x in (h_old, h_now, aice, area))
    if not (h_old.shape == h_now.shape == aice.shape == area.shape):
        raise ValueError(f"shape mismatch: h_old {h_old.shape}, h_now {h_now.shape}, "
                         f"aice {aice.shape}, area {area.shape}")
    if not dt_years > 0:
        raise ValueError(f"dt_years must be positive, got {dt_years!r}")
    area = np.where(np.isfinite(area) & (area > 0), area, 0.0)
    if area.sum() <= 0:
        raise ValueError("no cell has positive area")
    growth = (h_now - h_old) / dt_years
    S = growth * 0.5 * (h_old + h_now)
    full = (area > 0) & (aice >= config.full_cover) & (h_now > 0)
    if not full.any():
        raise ValueError(f"no fully ice-covered cell (aice >= {config.full_cover})")
    thick = full & (h_now >= np.median(h_now[full]))
    S_ref = float(np.median(S[thick]))
    if not S_ref > 0:
        raise ValueError(f"reference Stefan product S_ref = {S_ref:.3g} m^2/yr: the "
                         f"thick ice is not growing — no conduction-limited ice to jump")
    w = np.where(full, np.clip(S / S_ref, 0.0, 1.0), 0.0)
    weight = np.clip((w - config.ramp_lo) / (config.ramp_hi - config.ramp_lo), 0.0, 1.0)
    weight = np.where(full, weight, 0.0)
    return TaperMask(weight, w, growth, h_now, area, S_ref, float(dt_years), config)


def cell_factors(peak: float, weight) -> np.ndarray:
    """Per-cell factor ``1 + (peak - 1) * weight``."""
    return 1.0 + (float(peak) - 1.0) * np.asarray(weight, dtype=float)


def imbalance_ratio(peak: float, mask: TaperMask) -> float:
    """R(F) = (N_after - a) / (N_now - a) for peak factor ``peak``."""
    A, g = mask.area, mask.growth
    den = float((A * g).sum())
    if not den > 0:
        raise ValueError("area-integrated pre-jump ice growth is not positive: "
                         "the imbalance is not being carried by freezing")
    return float((A * g / cell_factors(peak, mask.weight)).sum() / den)


def effective_factor(peak: float, mask: TaperMask) -> float:
    """Ratio of area-mean ice after / before — what global-mean hi does."""
    A, h = mask.area, mask.h_now
    return float((A * h * cell_factors(peak, mask.weight)).sum() / (A * h).sum())


def solve_peak(r_target: float, mask: TaperMask, max_factor: float,
               search_max: float = 100.0) -> Optional[float]:
    """Peak factor F with R(F) = r_target, or None if unreachable by ``search_max``.

    R decreases monotonically from 1 at F = 1. ``max_factor`` is not applied
    here; the caller clips (and reports the unclipped value).
    """
    if not (0.0 < r_target < 1.0):
        raise ValueError(f"target ratio must be in (0, 1), got {r_target!r}")
    lo, hi = 1.0, float(search_max)
    if imbalance_ratio(hi, mask) > r_target:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if imbalance_ratio(mid, mask) > r_target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    return 0.5 * (lo + hi)
