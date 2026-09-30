"""aqua_ice plugin: scale the sea-ice reservoir in cice.r by one factor.

Replaces the prototype's ``accelerate.py --aqua_ice`` guess-and-check scaling
with a factor sized by the advisor and a constraint pass. Pure arrays; the
netCDF read/write lives in ``restart``.

What is scaled, and why this is self-consistent:

* ``vicen`` (ice volume per category) and ``eicen`` (ice enthalpy per layer)
  are multiplied by the same factor, so enthalpy per unit volume — the ice
  temperature/brine state — is unchanged; the ice simply gets thicker.
* ``aicen`` (ice area per category) is NOT touched. Ice fraction, surface
  temperature (``Tsfcn``) and hence albedo are unchanged, so the coupler
  exchange snapshot in cpl.r / cam.rs stays consistent with the ice state.
  This is what makes an in-place continuation restart safe for this plugin
  (restart-integration doc, O2).
* The ice factor may be one number or a per-cell ``(nj, ni)`` map (a tapered
  jump, ``taper.py``); a map multiplies every category and layer of a cell
  alike, so each cell keeps its own enthalpy per unit volume.
* Snow (``vsnon``/``esnon``) is scaled only when a separate snow factor != 1
  is requested; the Tier-0 sweep found snow near equilibrium already.
* Thickness per category rises (vicen/aicen); cells pushed past a category's
  upper bound are re-binned by CICE's own ITD remapping on the first step.

Water mass: in a slab-ocean aquaplanet the ocean is an unlimited water source,
so the added ice is not balanced against a finite inventory (unlike
Wordsworth's fixed-inventory case). The jump deliberately adds the ice the
skipped years would have frozen; the energy for it is the TOA deficit those
years would have radiated away.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from .plugins import PLUGIN_REGISTRY, ConstraintReport, VariablePlugin

ICE_FIELDS = ("vicen", "eicen")
SNOW_FIELDS = ("vsnon", "esnon")

#: Hard bound on any single factor, whatever the advisor says.
MAX_FACTOR = 2.0


def check_factor(factor: float, what: str = "factor") -> float:
    factor = float(factor)
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError(f"{what} must be a positive finite number, got {factor!r}")
    if not (1.0 / MAX_FACTOR <= factor <= MAX_FACTOR):
        raise ValueError(f"{what} {factor:g} outside the hard bound "
                         f"[{1.0 / MAX_FACTOR:g}, {MAX_FACTOR:g}]")
    return factor


def check_factor_map(factors, what: str = "factor") -> np.ndarray:
    """A per-cell factor map: every value within the hard bound."""
    f = np.asarray(factors, dtype=float)
    if f.ndim != 2:
        raise ValueError(f"{what} map must be 2-D (nj, ni), got shape {f.shape}")
    if not np.all(np.isfinite(f)):
        raise ValueError(f"{what} map has non-finite values")
    check_factor(f.min(), f"{what} (min)")
    check_factor(f.max(), f"{what} (max)")
    return f


class AquaIcePlugin(VariablePlugin):
    name = "aqua_ice"
    file_kind = "cice.r"
    variables = ICE_FIELDS + SNOW_FIELDS

    def apply_delta(self, fields: Dict[str, np.ndarray], delta) -> Dict[str, np.ndarray]:
        """``delta`` is ``(ice_factor, snow_factor)`` or a bare ice factor.

        ``ice_factor`` may be a ``(nj, ni)`` map; it broadcasts over the
        leading category/layer axis of ``vicen``/``eicen``.

        Only the fields present in ``fields`` are returned; snow fields are
        returned (scaled) only when the snow factor differs from 1.
        """
        if isinstance(delta, (tuple, list)):
            ice_f, snow_f = delta
        else:
            ice_f, snow_f = delta, 1.0
        if np.ndim(ice_f) == 0:
            ice_f = check_factor(ice_f, "ice factor")
        else:
            ice_f = check_factor_map(ice_f, "ice factor")
        snow_f = check_factor(snow_f, "snow factor")
        out = {}
        for name in ICE_FIELDS:
            out[name] = np.asarray(fields[name], dtype=float) * ice_f
        if snow_f != 1.0:
            for name in SNOW_FIELDS:
                out[name] = np.asarray(fields[name], dtype=float) * snow_f
        return out

    def enforce_constraints(
        self,
        fields_before: Dict[str, np.ndarray],
        fields_after: Dict[str, np.ndarray],
    ) -> Tuple[Dict[str, np.ndarray], ConstraintReport]:
        """Sign constraints, no-ice-from-nothing, and a totals report.

        Volumes must be >= 0 and enthalpies <= 0 (CICE stores energy of
        melting as negative). Scaling a valid state by a positive factor
        preserves both, so any violation here means the input was already
        unphysical — it is clamped and reported, not silently passed on.
        """
        out = dict(fields_after)
        adj = {}
        for name, arr in fields_after.items():
            arr = np.array(arr, dtype=float)
            before = np.asarray(fields_before[name], dtype=float)
            if name.startswith("v"):
                bad = arr < 0
                arr[bad] = 0.0
            else:
                bad = arr > 0
                arr[bad] = 0.0
            created = (before == 0) & (arr != 0)
            arr[created] = 0.0
            out[name] = arr
            tb, ta = float(before.sum()), float(arr.sum())
            ratio = ta / tb if tb != 0 else float("nan")
            msg = f"total {tb:.6g} -> {ta:.6g} (x{ratio:.4f})"
            if bad.any():
                msg += f"; clamped {int(bad.sum())} sign-violating cells"
            if created.any():
                msg += f"; zeroed {int(created.sum())} cells created from zero"
            adj[name] = msg
        return out, ConstraintReport(plugin=self.name, adjustments=adj)


PLUGIN_REGISTRY[AquaIcePlugin.name] = AquaIcePlugin()
