# CLAUDE.md — exocam-accelerate

Guidance for Claude Code sessions in this repository.

## What this is

A scientific tool to accelerate convergence of ExoCAM/CESM climate simulations
by extrapolating slow-evolving restart-file fields toward equilibrium
(`X_new = X + <dX/dt> * Δt`), generalizing Wordsworth et al. 2013 (Icarus,
§2.3) and Turbet et al. 2021 (Nature, Methods §4). It is a safe redesign of a
naive prototype (`accelerate.py`, kept at repo root for reference — it applies
a bare multiplicative factor to restart fields with none of the safeguards).

## Settled vs. open — respect this boundary

### Settled (implemented in `src/exocam_accelerate/`)

- Forward-Euler extrapolation of a measured tendency over a multi-year Δt.
- Tendencies come from **horizontally averaged per-layer trends**, not
  pointwise gridcell trends (pointwise amplifies dynamical noise). The fit
  helpers are shape-generic over trailing axes, so a plugin *may* fit
  pointwise trends (Wordsworth's surface-ice case is pointwise), but the
  default pipeline is per-layer means.
- **Trustworthiness gate before any step**: refuse — do not shrink — the step
  when trend-window curvature is large or the tendency sign is unstable.
  Large curvature means proximity to a nonlinear feedback threshold
  (ice-albedo, runaway greenhouse) where linear extrapolation is invalid.
- **Hard per-step magnitude clip**, configurable per variable (Turbet: 50 K).
- **Decreasing Δt schedule** (Wordsworth: 5×100 yr then 15×10 yr).
- Post-jump physical constraints are **per-variable plugins** (interface in
  `plugins.py`): e.g. aqua-ice enforces ice ≥ 0 / enthalpy ≤ 0 and no ice
  from nothing. (Wordsworth's water-mass conservation does not apply to a
  slab-ocean aquaplanet — the ocean is an unlimited source; see
  `aqua_ice.py`.) This is the only variable-specific piece.
- Post-step consistency check is a **hook** (interface in `consistency.py`):
  after the model runs forward from an accelerated restart, its own tendency
  must have continued along the predicted trajectory; otherwise roll back.
- **Cold-aquaplanet jump (decided 2026-09-25):** in-place edit of `cice.r`
  (`vicen`+`eicen`, one factor; `aicen` untouched) with a pristine
  `.pre-accel.nc` backup, then a plain continuation. Jump size comes from the
  Stefan conduction law `N = a + b/hi` (phase-space, target a chosen N, not
  N=0) — `advise.py`, `aqua_ice.py`, `restart.py`, `cli.py`. In-flight
  production safety (2026-09-26): pre-flight, `check`, `rollback` from the
  archived restart set (`runstate.py`, `check.py`); bookkeeping lives in
  `run/exocam_accelerate/` because CESM's st_archive sweeps
  `${CASE}.cice.r.*` in the run directory. See
  `docs/restart-integration-questions.md` §7 and
  `docs/phase-space-extrapolation.md`.

### Open (do NOT implement without a user decision)

- **Other restart paths/targets.** Hybrid reboot (`.i.` files), `somtp`
  (docn.r) and clm plugins, atmosphere fields. Remaining questions (O3–O8)
  are in `docs/restart-integration-questions.md`.
- Run orchestration (submitting/monitoring the forward runs between jumps).
- Turbet-style in-situ radiative heating-rate multiplication
  (`(P/Plim)^α`, α: 0.5 → 0.3 → 0) is **deferred indefinitely** — it requires
  modifying GCM Fortran source, out of scope for this file-manipulation tool.
  Recorded so the design space is not forgotten.

## Dependency boundary

- **`exocam-trend`** (`../exocam-trend`) is a **dependency** — but it is a flat
  script repo with no packaging, so we do not `import` from it. The stable
  interface is its `data/*.txt` time-series output (whitespace-separated
  columns: `month  VAR_native  VAR_int1  VAR_int2 ...`), parsed by
  `trend_io.py`. Its `trend_core.build_area_weights` was adapted (with
  attribution) into `trends.py` because we need per-layer horizontal means of
  3D fields, which exocam-trend does not produce (it only does global means of
  2D fields). If exocam-trend gains packaging or per-layer output, revisit.
- **`exocam-casemgr`** (`../exocam-casemgr`) is **reference-only**: read its
  restart handling and branch/hybrid setup to understand mechanics, but never
  import from it or couple to its internals. This tool will later be wired
  into casemgr as an external; keep the boundary clean.

## Layout

```
src/exocam_accelerate/
  trends.py       area weights, per-layer horizontal means, tendency/curvature fits
  safeguards.py   trustworthiness gate (refuse), hard magnitude clip
  schedule.py     decreasing Δt schedule
  stepper.py      gate → fit → extrapolate → clip pipeline (pure arrays)
  trend_io.py     parser for exocam-trend data/*.txt output (only file I/O in core)
  hindcast.py     Tier-0 offline hindcast harness
  phase_space.py  X-vs-N fits (linear / saturating / hyperbolic) + phase-space gate
  advise.py       jump advisor: trend columns -> ice factor, N_after, TS reference
  plugins.py      VariablePlugin interface + registry
  aqua_ice.py     aqua_ice plugin (scale vicen/eicen, constraint pass)
  restart.py      netCDF layer: cice.r in-place jump, backup/log in run/exocam_accelerate/
  runstate.py     pre-flight (SLURM, rpointers, archive) and whole-set rollback
  check.py        post-jump verdict PASS/WAIT/FAIL (the consistency check for aqua_ice)
  cli.py          `exocam-accelerate advise|jump|check|rollback|restore`
  consistency.py  generic ConsistencyCheck interface (design only; check.py is the aqua_ice one)
tests/            pytest unit tests for everything implemented
docs/restart-integration-questions.md   task-4 findings + open user decisions
```

## Conventions

- Time unit is **years** throughout; tendencies are per-year; Δt in years.
- Pure-computation modules must stay free of file I/O and netCDF imports;
  `netCDF4` belongs only to `restart.py` (optional extra `[netcdf]`).
- Reference restart files under `restart_set_examples/` are hundreds of MB:
  **header/metadata reads only** (`ncdump -h`), never read data payloads
  unless the user explicitly asks.
- `pytest` from the repo root runs the suite; keep it green.
