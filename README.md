# exocam-accelerate

Safeguarded convergence acceleration for ExoCAM/CESM climate simulations.

Some planetary climate regimes — snowball states with slowly growing ice
sheets, cold tidally locked planets, thick steam atmospheres — take hundreds
to thousands of simulated years to reach top-of-atmosphere equilibrium, at
prohibitive computational cost. This tool accelerates convergence by
extrapolating slow-evolving restart-file fields toward equilibrium using
tendencies measured from the simulation's own time series:

```
X_new = X + <dX/dt> * Δt
```

where `<dX/dt>` is a measured trend and `Δt` is a multi-year acceleration
timestep. The simulation is then restarted from the extrapolated state and run
forward, and the cycle repeats with a decreasing `Δt` until equilibrium.

## Scientific precedents

The method generalizes two published forward-Euler acceleration schemes:

- **Wordsworth, R., Forget, F., Millour, E., Head, J. W., Madeleine, J.-B., &
  Charnay, B. (2013)**, "Global modelling of the early martian climate under a
  denser CO2 atmosphere: Water cycle and ice evolution", *Icarus* 222, 1–19,
  §2.3 — ice-sheet equilibration by extrapolating the annual-mean ice rate of
  change with a 100 yr → 10 yr timestep schedule, enforcing ice ≥ 0 and total
  water-mass conservation after each jump.

- **Turbet, M., Bolmont, E., Chaverot, G., Ehrenreich, D., Leconte, J., &
  Marcq, E. (2021)**, "Day–night cloud asymmetry prevents early oceans on
  Venus but not on Earth", *Nature* 598, 276–280, Methods §4 — deep-atmosphere
  temperature convergence by extrapolating the **horizontally averaged**
  per-layer temperature tendency (`T_ijk += n_days × ΔT_mean,k`) with an
  n = 50 → 10 → 0 schedule and a hard 50 K per-step limit.

## What this package provides

The reusable core is the set of **safeguards** around the extrapolation, not
the per-variable scaling. Implemented as pure functions on arrays and trend
structures (no file I/O, no CESM coupling):

1. **Per-layer horizontal-mean tendencies** (`trends`) — tendencies are
   measured from horizontally averaged per-layer trends, not pointwise
   gridcell trends, isolating the slow coherent drift toward equilibrium from
   dynamical noise.
2. **Trustworthiness gate** (`safeguards`) — before any step, the trend window
   must show small curvature (d²X/dt²) and a stable tendency sign. Large
   curvature signals proximity to a nonlinear climate-feedback threshold
   (ice-albedo, runaway greenhouse) where linear extrapolation is unsafe; the
   gate **refuses** the step rather than shrinking it.
3. **Hard per-step magnitude clip** (`safeguards`) — an absolute per-variable
   limit on the step size, independent of what the extrapolation suggests
   (cf. Turbet's 50 K clip).
4. **Decreasing Δt schedule** (`schedule`) — large jumps early, small jumps
   near equilibrium, as in both precedents.

Two further safeguards are **interface-designed but not yet implemented**:

5. **Post-jump physical-constraint plugins** (`plugins`) — per-variable
   re-imposition of physical constraints after each jump (e.g. ice ≥ 0, total
   water-mass conservation). The only variable-specific piece of the design.
6. **Post-step consistency check** (`consistency`) — after the model runs
   forward from an accelerated restart, verify its own tendency continued
   along the predicted trajectory; if it corrected away, the step crossed a
   nonlinearity and must be rolled back.

The CESM restart-file write/run-continuation layer is **not implemented yet**;
see `docs/restart-integration-questions.md` for the open design questions.

## Related tools

- [`exocam-trend`](../exocam-trend) (dependency) — produces global-mean
  time-series diagnostics from monthly model output; its text output is one
  input pathway for tendency measurement (`trend_io`).
- [`exocam-casemgr`](../exocam-casemgr) (reference only) — case build/run/data
  management. This tool will eventually be invoked from casemgr as an
  external, but does not import from it.

## Install (development)

```bash
pip install -e ".[dev]"
pytest
```

## Status

Early scaffolding. Pure-computation safeguards implemented and tested;
restart-file manipulation and run orchestration intentionally absent.

**Use with care and caution — model behavior is not always straightforward.**
