"""Ocean jump advisor: from exocam-trend output to a somtp increment.

The hot / thick-atmosphere counterpart of ``advise`` (ice). Pure computation
over exocam-trend columns (read via ``trend_io``).

Regime: ice-free slab-ocean runs whose slow drift is the surface temperature
(docn.r ``somtp``) and the atmosphere coupled to it. There is no reservoir
with a law of its own like the Stefan ice: the whole surface-atmosphere system
relaxes toward N = 0 along its Gregory line, so the phase-space picture is
the classic one:

* ``TS = c0 + c1*N`` over a fit window (linear; the saturating form is fitted
  alongside, and a disagreement between the two is refused as a kink —
  near-runaway curvature). ``c1 < 0`` is required: TS must rise as N falls,
  otherwise there is no stable equilibrium in reach.
* The jump removes a fraction of the current imbalance: ``N_target =
  N_now*(1 - n_fraction)`` (or an explicit ``N_target`` between N_now and 0),
  so ``dTS = c1*(N_target - N_now)``.

Which imbalance. The slab ocean is in equilibrium when the net surface flux
into it vanishes (``energy_bot`` = 0; q-flux integrates to zero), and in the
hot runs it is the only reservoir the trend data can see: ``energy_bot`` =
C dTS/dt with C = the 50 m mixed layer and no offset (atlasfu D1-D5, 2026-09-30).
``energy_top`` differs from it by a near-constant offset that does not track
the warming rate — an atmospheric energy leak, up to ~10 W/m2 in the 4-bar
runs near 370 K — so ``energy_top`` -> 0 is not the equilibrium and a Gregory
line in it would overshoot by 5-15 K. The default coordinate is therefore
``energy_bot``; the leak (``energy_top - energy_bot`` over the window) is
reported.

Where the heat goes. Only somtp is edited. Whatever heat the atmosphere must
store to follow the surface (dry enthalpy plus, at 340-370 K, the vapor
column) is taken back from the mixed layer in the first weeks, so the run
lands short of the target by C_ocean / (C_ocean + C_atm). The trend data do
not resolve C_atm (it is collinear with the leak's temperature dependence),
so the default ``heat_ratio`` = 1 jumps somtp by the TS change; a first jump's
landing (``check``) measures the ratio for the next one.

Reported for the post-jump check and for judging the value of a jump:

* ``N_after``, ``TS_after`` — where the run should sit on the Gregory line;
* ``C_eff = N / (dTS/dt)`` over the window — the effective heat capacity the
  run is showing (also in metres of sea water, to compare with ``hblt``);
* the leak, mean ``energy_top - energy_bot`` over the window;
* ``years_skipped = tau*ln(N_now/N_after)`` with ``tau = -c1*C_eff`` — the
  one-box relaxation time the run would have needed.

After a jump (``since_year``, normally read from the jump log by the CLI):
only native annual means from ``settle_years`` after the jump are used, as for
the ice advisor. Automatic detection from the TS series is available but off
by default (see ``OceanAdvisorConfig.detect_jumps``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .advise import _annual, _have, model_years
from .phase_space import PhaseExtrapolation, PhaseGateConfig, extrapolate, fit_linear

PLUGIN = "som_ocean"
SURFACE_VAR = "TS"
SURFACE_FLUX_VAR = "energy_bot"
ICE_AREA_VAR = "ICEFRAC"

#: ``OceanAdvice.to_dict()`` shape. ``runstate.preflight`` refuses unknown ones.
OCEAN_ADVICE_SCHEMA_VERSION = "ocean-1"
#: advice refined by ``pattern_advice`` (per-cell increments; adds ``pattern``)
OCEAN_PATTERN_SCHEMA_VERSION = "ocean-2"

#: seconds per year and rho*cp of sea water (shr_const: 1026 kg/m3, 3996 J/kg/K)
_SEC_PER_YR = 365.0 * 86400.0
_RHO_CP_SW = 1026.0 * 3996.0


@dataclass(frozen=True)
class OceanAdvisorConfig:
    window_years: Optional[float] = None
    #: tried longest first; the longest whose TS(N) relation passes the gate
    auto_windows: Tuple[float, ...] = (40.0, 30.0, 20.0, 10.0)
    which: str = "int2"
    imbalance: str = "energy_bot"
    n_fraction: float = 0.5
    N_target: Optional[float] = None
    #: hard clip on the area-mean somtp increment, K
    max_dT: float = 10.0
    #: below max(min_dT, noise_factor x interannual TS scatter) a jump is not
    #: worth a restart edit (it would be lost in the year-to-year swings), K
    min_dT: float = 0.2
    noise_factor: float = 1.0
    #: somtp increment / TS change: 1 = no allowance for the heat the
    #: atmosphere takes back from the ocean after the jump
    heat_ratio: float = 1.0
    #: the regime is ice-free: refuse above this ICEFRAC anywhere in the window
    max_icefrac: float = 0.001
    #: the current state (mean of the last ``recent_years`` native years) must
    #: lie on the fitted line: |N_recent - N_line(TS_recent)| <= max(this,
    #: 2 sigma/sqrt(recent_years)). A steepening TS(N) (rising sensitivity,
    #: atlasfu D4 2026-09-30) otherwise puts a long window's line past the
    #: current state and advises a jump the wrong way.
    max_offline: float = 0.5
    recent_years: int = 5
    #: first model years never used (fast initial adjustment)
    spinup_years: int = 5
    gate: PhaseGateConfig = field(default_factory=PhaseGateConfig)

    since_year: Optional[int] = None
    #: off by default: hot runs swing 1-2 K a year and decelerate hard early,
    #: which a TS-step detector mistakes for jumps. The CLI takes the jump year
    #: from the jump log instead (``advise-ocean --rundir``). Forgetting it is
    #: benign: a jump moves the run along its Gregory line, and averaging
    #: points on a line keeps them on it.
    detect_jumps: bool = False
    #: a one-year |dTS| above this, standing out from its neighbours, is a jump
    jump_detect_dT: float = 1.0
    jump_detect_k: float = 4.0
    settle_years: int = 2
    post_jump_min_years: int = 15
    post_jump_gate: PhaseGateConfig = field(
        default_factory=lambda: PhaseGateConfig(min_abs_corr=0.7))

    def __post_init__(self) -> None:
        if not (0.0 < self.n_fraction <= 1.0):
            raise ValueError(f"n_fraction must satisfy 0 < n_fraction <= 1, got "
                             f"{self.n_fraction!r}")
        if not (0.0 < self.max_dT <= 25.0):
            raise ValueError(f"max_dT must satisfy 0 < max_dT <= 25 K, got "
                             f"{self.max_dT!r}")
        if self.min_dT < 0:
            raise ValueError("min_dT must be non-negative")
        if not (1.0 <= self.heat_ratio <= 3.0):
            raise ValueError(f"heat_ratio must satisfy 1 <= heat_ratio <= 3, got "
                             f"{self.heat_ratio!r}")
        if self.imbalance not in ("energy_bot", "energy_top"):
            raise ValueError(f"imbalance must be energy_bot or energy_top, got "
                             f"{self.imbalance!r}")
        if self.window_years is not None and self.window_years <= 0:
            raise ValueError(f"window_years must be positive, got {self.window_years!r}")
        if any(w <= 0 for w in self.auto_windows):
            raise ValueError(f"auto_windows must be positive, got {self.auto_windows!r}")
        for name in ("settle_years", "post_jump_min_years"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.spinup_years < 0:
            raise ValueError("spinup_years must be non-negative")


@dataclass(frozen=True)
class OceanAdvice:
    case: str
    model_year: int
    N_now: float                  # native annual mean, last year
    N_now_fit: float              # Gregory line at TS_now
    N_target: float
    TS_now: float
    config: OceanAdvisorConfig
    which_used: str
    window_years: float
    since_year: Optional[int]
    detected_jumps: Tuple[int, ...]
    gregory: Optional[PhaseExtrapolation]
    leak: dict                    # mean/trend of energy_top - energy_bot, window
    dTS_target: Optional[float]   # TS change the target asks for
    somtp_dT: Optional[float]     # None -> do not jump
    somtp_dT_raw: Optional[float]
    clipped: bool
    TS_after: Optional[float]
    N_after: Optional[float]
    C_eff: Optional[float]        # W yr m-2 K-1
    tau_years: Optional[float]
    years_skipped: Optional[float]
    reasons: tuple
    warnings: tuple
    provenance: Optional[Dict[str, str]] = None
    #: every window tried: {window, c0, c1, corr, offline, status}
    windows_tried: tuple = ()

    @property
    def jump(self) -> bool:
        return self.somtp_dT is not None

    def to_dict(self) -> dict:
        def num(x):
            return None if x is None or not np.isfinite(x) else float(x)

        g = self.gregory
        greg = None
        if g is not None:
            greg = {"form": g.form, "accepted": g.accepted, "reasons": list(g.reasons),
                    "N_target": num(g.N_target), "X_now": num(g.X_now),
                    "prediction": num(g.prediction),
                    "predictions": {k: num(v) for k, v in g.predictions.items()},
                    "fits": {k: {"params": list(f.params), "corr": f.corr, "rmse": f.rmse}
                             for k, f in g.fits.items()},
                    "extrapolation_ratio": num(g.extrapolation_ratio)}
        c = self.config
        C_m = None if self.C_eff is None else self.C_eff * _SEC_PER_YR / _RHO_CP_SW
        return {
            "schema_version": OCEAN_ADVICE_SCHEMA_VERSION,
            "plugin": PLUGIN,
            "case": self.case,
            "model_year": self.model_year,
            "N_now": num(self.N_now),
            "N_now_fit": num(self.N_now_fit),
            "N_target": num(self.N_target),
            "TS_now": num(self.TS_now),
            "which": self.which_used,
            "window_years": self.window_years,
            "since_year": self.since_year,
            "detected_jumps": list(self.detected_jumps),
            "config": {"n_fraction": c.n_fraction, "N_target": c.N_target,
                       "max_dT": c.max_dT, "min_dT": c.min_dT,
                       "noise_factor": c.noise_factor, "max_offline": c.max_offline,
                       "heat_ratio": c.heat_ratio,
                       "imbalance": c.imbalance, "max_icefrac": c.max_icefrac},
            "gregory": greg,
            "leak": {k: num(v) for k, v in self.leak.items()},
            "dTS_target": num(self.dTS_target),
            "somtp_dT": num(self.somtp_dT),
            "somtp_dT_raw": num(self.somtp_dT_raw),
            "somtp_dT_clipped": self.clipped,
            "TS_after": num(self.TS_after),
            "N_after": num(self.N_after),
            "C_eff_W_yr_m2_K": num(self.C_eff),
            "C_eff_m_seawater": num(C_m),
            "tau_years": num(self.tau_years),
            "years_skipped": num(self.years_skipped),
            "expected_after_jump": {"TS": num(self.TS_after)},
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "provenance": dict(self.provenance) if self.provenance else None,
            "windows_tried": [{k: (num(v) if isinstance(v, float) else v)
                               for k, v in d.items()} for d in self.windows_tried],
        }


def detect_ts_jumps(years, ts, dT: float = 1.0, k_noise: float = 4.0,
                    before: int = 6, after: int = 2) -> List[int]:
    """Model years at which annual-mean TS steps off its own recent trajectory.

    For each candidate year, a quadratic is fitted to the ``before`` preceding
    years and extrapolated over the next ``after``; the candidate is a jump
    when every one of those years departs, with one sign, by more than both
    ``dT`` and ``k_noise`` x the residual scatter of that fit. Hot runs swing by 1-2 K from year to year and warm by several
    K/yr early on, so a bare one-year step test fires on both; a jump is a
    persistent offset from the trend. Returns the first model year run from
    each jumped state.
    """
    ts = np.asarray(ts, dtype=float)
    years = np.asarray(years)
    out: List[int] = []
    for i in range(before, ts.size - after + 1):
        x = np.arange(i - before, i)
        y = ts[i - before:i]
        if not np.all(np.isfinite(y)) or not np.all(np.isfinite(ts[i:i + after])):
            continue
        c = np.polyfit(x, y, 2)
        sig = float(np.std(y - np.polyval(c, x), ddof=3))
        deps = ts[i:i + after] - np.polyval(c, np.arange(i, i + after))
        thr = max(dT, k_noise * sig)
        if np.all(np.abs(deps) > thr) and abs(np.sign(deps).sum()) == deps.size:
            if out and years[i] - out[-1] <= after:      # one jump, one report
                continue
            out.append(int(years[i]))
    return out


def _leak(columns, which, mask) -> dict:
    """``energy_top - energy_bot`` over the window: mean (W/m2) and trend (W/m2/yr)."""
    if not (_have(columns, "energy_top", which) and _have(columns, "energy_bot", which)):
        return {"mean": None, "trend_per_yr": None}
    t, top = _annual(columns, "energy_top", which)
    _, bot = _annual(columns, "energy_bot", which)
    g = (top - bot)[mask]
    ok = np.isfinite(g)
    if ok.sum() < 2:
        return {"mean": None, "trend_per_yr": None}
    return {"mean": float(g[ok].mean()),
            "trend_per_yr": float(np.polyfit(t[mask][ok], g[ok], 1)[0])}


def advise_ocean(columns: Dict[str, np.ndarray], case: str,
                 config: OceanAdvisorConfig = OceanAdvisorConfig(),
                 start_year: int = 1,
                 provenance: Optional[Dict[str, str]] = None) -> OceanAdvice:
    """Build a somtp jump advice from merged exocam-trend columns (``load_case``)."""
    t, N_nat = _annual(columns, config.imbalance, "native")
    years = model_years(t, start_year)
    model_year = int(years[-1])
    N_now = float(N_nat[-1])
    reasons: List[str] = []
    warnings: List[str] = []
    windows: List[dict] = []
    have_T = _have(columns, SURFACE_VAR, "native")

    detected: List[int] = []
    TS_nat = None
    if have_T:
        _, TS_nat = _annual(columns, SURFACE_VAR, "native")
        if config.detect_jumps:
            detected = detect_ts_jumps(years, TS_nat, config.jump_detect_dT,
                                       config.jump_detect_k)
    since = config.since_year
    if since is None and detected:
        since = detected[-1]
        warnings.append(f"detected a jump in {SURFACE_VAR} at model year {since}; "
                        f"using post-jump data only (override with --since)")
    elif since is not None and detected and detected[-1] > since:
        warnings.append(f"a later jump is visible in {SURFACE_VAR} at model year "
                        f"{detected[-1]} than --since {since}")

    if since is not None:
        which, gate = "native", config.post_jump_gate
        usable = years >= since + config.settle_years
    else:
        which, gate = config.which, config.gate
        usable = years >= years[0] + config.spinup_years
    # fit only the ice-free part of the run (early years melt the initial ice)
    icy = None
    if _have(columns, ICE_AREA_VAR, "native"):
        _, ice = _annual(columns, ICE_AREA_VAR, "native")
        icy = ~(ice <= config.max_icefrac)            # NaN counts as icy
        if icy.any():
            usable &= years > years[icy].max()
    cap = float(usable.sum())

    def result(N_now_fit=float("nan"), N_target=float("nan"), window=0.0, greg=None,
               leak=None, dTS=None, dT=None, raw=None, clipped=False, TS_after=None,
               N_after=None, C_eff=None, tau=None, skipped=None):
        return OceanAdvice(case, model_year, N_now, float(N_now_fit), float(N_target),
                           float(TS_nat[-1]) if TS_nat is not None else float("nan"),
                           config, which, float(window), since, tuple(detected), greg,
                           leak or {}, dTS, dT, raw, clipped, TS_after, N_after, C_eff,
                           tau, skipped, tuple(reasons), tuple(warnings), provenance,
                           tuple(windows))

    if not have_T:
        reasons.append(f"no {SURFACE_VAR} series in the trend output")
        return result()
    if since is not None and cap < config.post_jump_min_years:
        reasons.append(f"only {int(cap)} settled years since the jump at model year "
                       f"{since} (need {config.post_jump_min_years}, after "
                       f"{config.settle_years} settling years): keep running")
        return result()
    if icy is None:
        reasons.append(f"no {ICE_AREA_VAR} series: cannot confirm the run is ice-free "
                       f"(the ocean jump never warms water under ice)")
        return result()
    if icy[-1]:
        reasons.append(f"{ICE_AREA_VAR} is {ice[-1]:.4f} in the latest year (> "
                       f"{config.max_icefrac:g}): not the ice-free regime — use the "
                       f"aqua_ice advisor")
        return result()
    TS_now = float(TS_nat[-1])
    if not (np.isfinite(TS_now) and np.isfinite(N_now)):
        reasons.append("the latest annual TS or N is not finite")
        return result()

    _, N_fit = _annual(columns, config.imbalance, which)
    _, T_fit = _annual(columns, SURFACE_VAR, which)

    # ---- window: longest candidate whose TS(N) relation passes the gate ----
    if config.window_years is not None:
        cands = [min(float(config.window_years), cap)]
    else:
        cands = sorted({W for W in config.auto_windows if W <= cap}, reverse=True)
        if since is not None:
            cands = sorted({min(cap, max(config.auto_windows))} |
                           {W for W in cands if W >= config.post_jump_min_years},
                           reverse=True)
        if not cands:
            cands = [cap]

    def fit_window(W):
        m = usable & (t > t[-1] - W)
        f = None
        Nw, Tw = N_fit[m], T_fit[m]
        ok = np.isfinite(Nw) & np.isfinite(Tw)
        if ok.sum() >= 2 and np.ptp(Nw[ok]) > 0:
            f = fit_linear(Nw[ok], Tw[ok])
        return m, f

    Nn = N_nat
    k = min(config.recent_years, int(np.isfinite(Nn).sum()))
    N_rec = float(np.nanmean(Nn[-k:]))
    T_rec = float(np.nanmean(TS_nat[-k:]))
    tail = Nn[-2 * k:]
    tail = tail[np.isfinite(tail)]
    sig = float(np.std(tail - np.polyval(np.polyfit(np.arange(tail.size), tail, 1),
                                         np.arange(tail.size)))) if tail.size > 2 else 0.0
    tol_off = max(config.max_offline, 2.0 * sig / np.sqrt(max(k, 1)))
    ttail = TS_nat[-2 * k:]
    ttail = ttail[np.isfinite(ttail)]
    sig_T = float(np.std(ttail - np.polyval(np.polyfit(np.arange(ttail.size), ttail, 1),
                                            np.arange(ttail.size)))) if ttail.size > 2 else 0.0
    min_move = max(config.min_dT, config.noise_factor * sig_T)

    chosen = None
    tried: List[str] = []
    for W in cands:
        m, f = fit_window(W)
        if f is None:
            continue
        c0, c1 = f.params
        rec = {"window": float(W), "c0": float(c0), "c1": float(c1),
               "corr": float(f.corr), "offline": None, "status": "ok"}
        windows.append(rec)
        if c1 >= 0:
            rec["status"] = f"dTS/dN {c1:+.2f} >= 0"
            tried.append(f"{W:g} yr: {rec['status']}")
            continue
        N_now_fit = (TS_now - c0) / c1
        off = N_rec - (T_rec - c0) / c1
        rec["offline"] = float(off)
        probe = extrapolate(SURFACE_VAR, N_fit[m], T_fit[m], 0.5 * N_now_fit, TS_now,
                            "linear", gate)
        if not probe.accepted:
            rec["status"] = "; ".join(probe.reasons)
            tried.append(f"{W:g} yr: {rec['status']}")
            continue
        if abs(off) > tol_off:
            rec["status"] = (f"current state {off:+.2f} W/m2 off the line "
                             f"(tolerance {tol_off:.2f})")
            tried.append(f"{W:g} yr: {rec['status']}")
            continue
        chosen = (W, m, f)
        break
    if chosen is None:
        reasons.append("no fit window gives a TS(N) line that is trustworthy and "
                       "passes through the current state (curving relation — e.g. "
                       "sensitivity rising with temperature — or too noisy): "
                       + " | ".join(tried) if tried else
                       "no window with a usable TS(N) relation")
        return result(window=cands[0],
                      leak=_leak(columns, which, usable & (t > t[-1] - cands[0])))
    window, w, lin = chosen
    c0, c1 = lin.params
    N_now_fit = (TS_now - c0) / c1

    if config.N_target is not None:
        N_target = float(config.N_target)
        lo, hi = sorted((N_now_fit, 0.0))
        if not (lo <= N_target <= hi) or N_target == N_now_fit:
            reasons.append(f"explicit N_target {N_target:+.2f} is not between the "
                           f"current imbalance {N_now_fit:+.2f} and 0: refusing (would "
                           f"not advance toward equilibrium)")
    else:
        N_target = N_now_fit * (1.0 - config.n_fraction)
    Nw = N_fit[w][np.isfinite(N_fit[w])]
    reach = gate.max_extrapolation_ratio * float(np.ptp(Nw))
    if N_target > Nw.max() + reach or N_target < Nw.min() - reach:
        wanted = N_target
        N_target = float(Nw.max() + reach if N_target > Nw.max() else Nw.min() - reach)
        warnings.append(f"target N {wanted:+.2f} limited to {N_target:+.2f} by the "
                        f"window's N-range ({gate.max_extrapolation_ratio:g}x "
                        f"{np.ptp(Nw):.2f} W/m2)")
    greg = extrapolate(SURFACE_VAR, N_fit[w], T_fit[w], N_target, TS_now, "linear", gate)
    reasons.extend(greg.reasons)

    leak = _leak(columns, which, w)
    if (config.imbalance == "energy_top" and leak["mean"] is not None
            and abs(leak["mean"]) > 1.0):
        warnings.append(f"energy_top - energy_bot = {leak['mean']:+.1f} W/m2 over the "
                        f"window: energy_top -> 0 is not the ocean's equilibrium; "
                        f"prefer the energy_bot coordinate")

    # ---- effective heat capacity and relaxation time ----
    # The heat taken up since the window start, sum(N dt), against TS: its
    # slope is C for any trajectory (no d/dt of a noisy series, no bias from
    # curvature), and summing integrates out the interannual noise in N.
    C_eff = tau = None
    _, Nn_w = _annual(columns, config.imbalance, "native")
    ok = np.isfinite(Nn_w[w]) & np.isfinite(TS_nat[w])
    if ok.sum() >= 3:
        Nk = Nn_w[w][ok]
        E = np.cumsum(Nk) - 0.5 * Nk          # heat to mid-year, W yr m-2
        C = float(np.polyfit(TS_nat[w][ok], E, 1)[0])
        if C > 0:
            C_eff = C
            tau = -c1 * C_eff

    dTS = dT = raw = TS_after = N_after = skipped = None
    clipped = False
    if not reasons:
        dTS = float(greg.prediction - TS_now)
        raw = dTS * config.heat_ratio
        dT = float(np.clip(raw, -config.max_dT, config.max_dT))
        clipped = dT != raw
        if clipped:
            warnings.append(f"somtp increment {raw:+.2f} K clipped to {dT:+.2f} K")
        if abs(dT / config.heat_ratio) < min_move:
            reasons.append(f"the jump would move TS by only "
                           f"{dT / config.heat_ratio:+.2f} K (< {min_move:.2f} K: the "
                           f"larger of {config.min_dT:g} and the interannual TS "
                           f"scatter): not worth a restart edit — near equilibrium")
            dT = raw = None
            clipped = False
        else:
            TS_after = TS_now + dT / config.heat_ratio
            N_after = float((TS_after - c0) / c1)
            if tau is not None and N_now_fit * N_after > 0 and abs(N_after) < abs(N_now_fit):
                skipped = float(tau * np.log(N_now_fit / N_after))
    return result(N_now_fit, N_target, window, greg, leak, dTS if dT is not None else None,
                  dT, raw, clipped, TS_after, N_after, C_eff, tau, skipped)


def pattern_advice(advice: dict, summary: dict) -> dict:
    """Mark an ocean advice (``OceanAdvice.to_dict()``) as a patterned jump.

    The area-mean increment is unchanged (the weights have area mean 1), so
    the expected TS/N after the jump are unchanged; the result records the
    pattern (``summary`` from ``OceanPattern.summary()`` plus the weight-file
    bookkeeping the CLI adds) and bumps the schema.
    """
    if advice.get("plugin") != PLUGIN:
        raise ValueError("not an ocean (som_ocean) advice")
    if advice.get("somtp_dT") is None:
        raise ValueError("advice says do not jump; nothing to pattern")
    if advice.get("pattern"):
        raise ValueError("advice already carries a pattern")
    out = dict(advice)
    out.update(schema_version=OCEAN_PATTERN_SCHEMA_VERSION, pattern=dict(summary))
    return out


# ---------------------------------------------------------------------------
# probe mode: a warm (or cool) perturbation to map the trajectory
# ---------------------------------------------------------------------------
#
# Where the Gregory line is refused because the recent relation is not
# constrained (atlasfu D4/D5, 2026-09-30: over the last 10-20 yr energy_bot
# moves by less than its noise while TS rises 3-6 K, so lambda is anywhere in
# 0-0.5 W/m2/K and the equilibrium anywhere above ~375 K), the jump is not
# sized to an equilibrium. It is a probe: step TS ahead by the run's own recent
# warming over ``probe_years`` (forward Euler behind the time-domain
# trustworthiness gate and the hard clip), then read the response. The probe
# supplies the lever arm the natural run lacks: a settled post-probe point
# several K away from the pre-probe point measures lambda directly, and the
# sign of energy_bot after it says which side of the equilibrium the run is on
# (``check.check_ocean_probe``).

OCEAN_PROBE_SCHEMA_VERSION = "ocean-probe-1"


@dataclass(frozen=True)
class ProbeConfig:
    #: the probe steps TS by (recent dTS/dt) x probe_years
    probe_years: float = 15.0
    #: explicit probe size, K: skips the trend sizing and its gate (the
    #: user's judgment, e.g. a run whose year-to-year swings fail the gate)
    probe_dT: Optional[float] = None
    #: years of native annual TS the recent trend is fitted over
    trend_years: int = 15
    max_dT: float = 10.0
    #: the probe must exceed this many interannual TS sigmas (so it is visible)
    min_sigmas: float = 3.0
    #: years averaged for the pre-probe point and in the lambda estimate
    recent_years: int = 5
    imbalance: str = "energy_bot"
    max_icefrac: float = 0.001
    spinup_years: int = 5
    since_year: Optional[int] = None
    settle_years: int = 2

    def __post_init__(self) -> None:
        if not (0.0 < self.probe_years <= 50.0):
            raise ValueError(f"probe_years must be in (0, 50], got {self.probe_years!r}")
        if self.trend_years < 5:
            raise ValueError("trend_years must be at least 5")
        if not (0.0 < self.max_dT <= 25.0):
            raise ValueError(f"max_dT must satisfy 0 < max_dT <= 25 K, got {self.max_dT!r}")
        if self.recent_years < 3:
            raise ValueError("recent_years must be at least 3")
        if self.imbalance not in ("energy_bot", "energy_top"):
            raise ValueError(f"imbalance must be energy_bot or energy_top, got "
                             f"{self.imbalance!r}")


def probe_ocean(columns: Dict[str, np.ndarray], case: str,
                config: ProbeConfig = ProbeConfig(), start_year: int = 1,
                provenance: Optional[Dict[str, str]] = None) -> dict:
    """Probe advice (a JSON-able dict, ``somtp_dT`` None = do not probe)."""
    from .safeguards import GateConfig
    from .stepper import propose_step
    from .trends import TrendSeries

    t, T = _annual(columns, SURFACE_VAR, "native")
    years = model_years(t, start_year)
    reasons: List[str] = []
    warnings: List[str] = []
    k = config.recent_years
    out = {"schema_version": OCEAN_PROBE_SCHEMA_VERSION, "plugin": PLUGIN,
           "mode": "probe", "case": case, "model_year": int(years[-1]),
           "since_year": config.since_year,
           "config": {"probe_years": config.probe_years, "probe_dT": config.probe_dT,
                      "trend_years": config.trend_years, "max_dT": config.max_dT,
                      "min_sigmas": config.min_sigmas, "recent_years": k,
                      "imbalance": config.imbalance},
           "provenance": dict(provenance) if provenance else None}

    usable = years >= (config.since_year + config.settle_years
                       if config.since_year is not None
                       else years[0] + config.spinup_years)
    if not _have(columns, ICE_AREA_VAR, "native"):
        reasons.append(f"no {ICE_AREA_VAR} series: cannot confirm the run is ice-free")
    else:
        _, ice = _annual(columns, ICE_AREA_VAR, "native")
        icy = ~(ice <= config.max_icefrac)
        if icy[-1]:
            reasons.append(f"{ICE_AREA_VAR} is {ice[-1]:.4f} in the latest year: not "
                           f"the ice-free regime")
        elif icy.any():
            usable &= years > years[icy].max()
    w = usable & (t > t[-1] - config.trend_years)
    _, N = _annual(columns, config.imbalance, "native")
    if not reasons and w.sum() < config.trend_years:
        reasons.append(f"only {int(w.sum())} usable years for the trend (need "
                       f"{config.trend_years})")
    if not reasons and not (np.all(np.isfinite(T[w])) and np.all(np.isfinite(N[w]))):
        reasons.append("non-finite TS or N in the trend window")

    if not reasons:
        step = propose_step(TrendSeries(t[w], T[w]), config.probe_years, config.max_dT,
                            GateConfig())
        tr = np.polyfit(t[w], T[w], 1)
        sig_T = float(np.std(T[w] - np.polyval(tr, t[w]), ddof=2))
        Nr = N[w][-k:]
        sig_N = float(np.std(N[w] - np.polyval(np.polyfit(t[w], N[w], 1), t[w]), ddof=2))
        out.update(TS_now=float(T[-1]), trend_K_per_yr=float(tr[0]), sigma_TS=sig_T,
                   N_now=float(Nr.mean()), sigma_N=sig_N)
        # two-period lambda over the window (context: why a probe is needed)
        a, b = N[w][-2 * k:-k], N[w][-k:]
        dTp = float(T[w][-k:].mean() - T[w][-2 * k:-k].mean())
        if a.size == k and dTp != 0:
            se = float(np.sqrt(a.var(ddof=1) / k + b.var(ddof=1) / k))
            out["lambda_recent"] = {"value": float(-(b.mean() - a.mean()) / dTp),
                                    "se": abs(se / dTp), "dTS": dTp}
        if config.probe_dT is not None:
            dT = float(np.clip(config.probe_dT, -config.max_dT, config.max_dT))
            out["somtp_dT_raw"] = float(config.probe_dT)
            out["somtp_dT_clipped"] = dT != config.probe_dT
            out["sizing"] = "explicit"
            if not step.accepted:
                warnings.append("the trend gate would have refused an automatic probe: "
                                + "; ".join(step.gate.reasons))
        elif not step.accepted:
            reasons.extend(step.gate.reasons)
        else:
            dT = float(step.delta)
            out["somtp_dT_raw"] = float(step.tendency * config.probe_years)
            out["somtp_dT_clipped"] = bool(step.clip.any_clipped)
            out["sizing"] = "trend"
        if not reasons:
            if abs(dT) < config.min_sigmas * sig_T:
                reasons.append(f"probe {dT:+.2f} K is within {config.min_sigmas:g} "
                               f"interannual TS sigmas ({sig_T:.2f} K): its response "
                               f"would not be readable — lengthen --probe-years")
            else:
                n = k
                se_dN = sig_N * np.sqrt(2.0 / n)
                out.update(somtp_dT=dT, TS_after=float(T[-1] + dT),
                           expected_after_jump={"TS": float(T[-1] + dT)},
                           # smallest lambda the check can tell from zero
                           lambda_detectable=float(2.0 * se_dN / abs(dT)),
                           years_skipped=(float(config.probe_years) if out["sizing"] == "trend"
                               else float(abs(dT / tr[0])) if tr[0] * dT > 0 else None))
                if out["somtp_dT_clipped"]:
                    warnings.append(f"probe {out['somtp_dT_raw']:+.2f} K clipped to "
                                    f"{dT:+.2f} K")
    out.setdefault("somtp_dT", None)
    out["reasons"], out["warnings"] = reasons, warnings
    return out
