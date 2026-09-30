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
- **Tapered jump (decided 2026-09-29):** after the first real jumps (grp4
  pt01/pt03 ×1.5 at 0141, both PASS) showed near-equilibrium thin ice (the
  substellar ring) melting back, the jump is applied per cell as
  `1 + (F-1)·weight`, with weight a ramp (default 0.05–0.40) of the Stefan
  product `dh/dt·h` relative to the thick ice, measured between two pristine
  archived restarts; partial-ice cells are never scaled. Peak `F` solves
  `(N_after-a)/(N_now-a) = Σ A g/f / Σ A g` (per-cell freezing ∝ 1/h).
  `taper.py` (pure), `advise.taper_advice` (schema 2: `effective_factor`,
  `law_offset_after`), `restart.build_taper_mask`, `cli taper`; `check` scores
  against law + offset with the effective factor.

- **Hot, ice-free ocean jump (decided 2026-09-30, branch `feature/ocean-jump`):**
  in-place edit of `docn.r` `somtp` (K, lon-fastest `gsize = ni*nj` on the
  docn domain named in `docn_ocn_in`), same backup/log/pre-flight/check/
  rollback machinery, plain continuation. Phase-space coordinate is
  **`energy_bot`** (the slab's own equilibrium): `energy_top` carries an
  atmospheric energy leak of up to ~10 W/m² in the 4-bar atlasfu runs. Gregory
  line TS(N_bot), current state must lie on it (refuses the steepening D4/D5),
  jump floor = interannual TS scatter, `--heat-ratio` (default 1) for heat the
  atmosphere takes back. Uniform by default; `pattern` = measured warming
  pattern (user asked for both). `som_ocean.py` (pure), `ocean_advise.py`,
  `restart.apply_ocean_jump`/`write_somtp_map`, `check.check_ocean_jump`,
  `cli advise-ocean|pattern|somtp-map`, viewer ocean layout + maps. See
  `docs/ocean-jump.md`.

### Open (do NOT implement without a user decision)

- **Other restart paths/targets.** Hybrid reboot (`.i.` files), clm plugins,
  atmosphere fields (cam.r `PT` = T_v/pkz is editable for a continuation, but
  the pbuf copies — `TCWAT` above all — must follow; see
  `docs/ocean-jump.md`). Remaining questions (O3–O8)
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
  taper.py        tapered jump: per-cell Stefan weights, peak-factor solve
  aqua_ice.py     aqua_ice plugin (scale vicen/eicen by a factor or per-cell map)
  som_ocean.py    som_ocean plugin (shift docn.r somtp; warming-pattern weights)
  ocean_advise.py ocean jump advisor: TS vs energy_bot Gregory line -> somtp increment
  restart.py      netCDF layer: cice.r / docn.r in-place jumps, backup/log in
                  run/exocam_accelerate/, somtp lat-lon maps
  runstate.py     pre-flight (SLURM, rpointers, archive) and whole-set rollback
  check.py        post-jump verdict PASS/WAIT/FAIL (aqua_ice and som_ocean checks)
  cli.py          `exocam-accelerate advise|taper|advise-ocean|pattern|somtp-map|
                  jump|check|rollback|restore|view`
  viewer.py       local interactive viewer (stdlib http.server + viewer_assets/index.html,
                  Plotly via CDN); draws only what advise/check compute
  consistency.py  generic ConsistencyCheck interface (design only; check.py is the aqua_ice one)
tests/            pytest unit tests for everything implemented
docs/restart-integration-questions.md   task-4 findings + open user decisions
docs/ocean-jump.md                       hot-regime (somtp) jump: physics findings + design
```

## Conventions

- Time unit is **years** throughout; tendencies are per-year; Δt in years.
- Pure-computation modules must stay free of file I/O and netCDF imports;
  `netCDF4` belongs only to `restart.py` (optional extra `[netcdf]`).
- Reference restart files under `restart_set_examples/` are hundreds of MB:
  **header/metadata reads only** (`ncdump -h`), never read data payloads
  unless the user explicitly asks.
- `pytest` from the repo root runs the suite; keep it green.
