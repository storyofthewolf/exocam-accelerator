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
  `plugins.py`): e.g. aqua-ice must enforce ice ≥ 0 and conserve total water
  mass. This is the only variable-specific piece.
- Post-step consistency check is a **hook** (interface in `consistency.py`):
  after the model runs forward from an accelerated restart, its own tendency
  must have continued along the predicted trajectory; otherwise roll back.

### Open (do NOT implement without a user decision)

- **The CESM restart-write / run-continuation layer.** The central open
  question: apply an accelerated state by modifying `.r.` files and doing a
  simple continuation, or by modifying/creating `.i.` files and rebooting as a
  hybrid run — and how each path treats perturbed fields (read verbatim vs.
  re-derived at init). Findings and remaining decisions are in
  `docs/restart-integration-questions.md`. Read that before touching this.
- Implementations of the plugin and consistency interfaces.
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
  plugins.py      VariableConstraint interface (design only — NotImplementedError)
  consistency.py  ConsistencyCheck interface (design only — NotImplementedError)
tests/            pytest unit tests for everything implemented
docs/restart-integration-questions.md   task-4 findings + open user decisions
```

## Conventions

- Time unit is **years** throughout; tendencies are per-year; Δt in years.
- Pure-computation modules must stay free of file I/O and netCDF imports;
  `netCDF4` belongs only to the future restart layer (optional extra).
- Reference restart files under `restart_set_examples/` are hundreds of MB:
  **header/metadata reads only** (`ncdump -h`), never read data payloads
  unless the user explicitly asks.
- `pytest` from the repo root runs the suite; keep it green.
