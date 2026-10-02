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

## Offline validation (Tier-0 hindcast)

Before any restart file is touched, the method can be validated against ground
truth the model archives already contain. The `hindcast` module scores
forward-Euler extrapolation on existing time series with no new model runs:
truncate an archived series at year *N*, fit a tendency over a trailing window,
extrapolate by Δt through the real gate-and-clip pipeline, and compare to what
the simulation actually did at *N* + Δt. Sweeping (*N*, window, Δt) over a set
of archived runs yields:

- **skill maps** — forward-Euler error per variable / regime / Δt
  (`error_by_dt`);
- **calibrated gate thresholds** — of the steps the gate accepted, how many the
  extrapolation actually got right, and how many refusals were unnecessary
  (`calibrate_gate`: `miss_rate`, `false_alarm_rate`);
- **a defensible Δt schedule** — the largest Δt whose accepted hindcasts stay
  within tolerance (`max_safe_dt`), rather than a hand-tuned one.

This is Tier 0 of a five-tier test plan (`docs/restart-integration-questions.md`
§6b) that progresses from offline hindcasts to a full twin-convergence
experiment. It consumes `exocam-trend` text output via `trend_io`; the batch
driver that generates that output lives in the `exocam-trend` repo.

## Jumping a cold aquaplanet case (in-flight production runs)

For cold/waterbelt aquaplanets whose slow drift is sea-ice growth. The TOA
deficit `N` is conducted through the ice and freezes onto it, so
`N ≈ a + b/hi` (Stefan-limited growth; `docs/phase-space-extrapolation.md`).
A jump picks a target imbalance and sets the ice to the thickness that law
says produces it, by scaling ice volume and enthalpy (`vicen`, `eicen`) in
`cice.r` by one factor. Ice area, surface temperature and albedo are left
alone, so the run resumes as a plain continuation.

**Tapered jump (recommended).** A uniform factor also thickens ice that was
already near its local equilibrium — on tidally locked cases, the thin ice
around the substellar point — and that extra ice melts back (grp4 pt01/pt03,
2026-09-29). `taper` weights each cell by how conduction-limited its growth
was over the previous `--baseline-years` (Stefan product `dh/dt·h` against the
thick ice's), scales cell by cell as `1 + (F-1)·weight`, and solves for the
peak factor `F` that still reaches the target imbalance. It reads two
archived `cice.r` files and a `cice.h` grid file, so it runs on the HPC; the
tapered advice (schema 2) and its `.taper.nc` weight map feed `jump` and
`check` unchanged.

The workflow is built for babysitting a few production runs at a time: every
step that touches a case happens at a segment boundary with no job queued,
every jump can be undone, and the post-jump verdict is one command.

**Stage 0 — shadow (all cases, no risk).** At each segment boundary run
`advise` and keep the JSON. Only cases whose advice is stable across segments
are candidates.

**Stage 1 — one canary** (ideally a branch clone of a production case), jumped
at a gentle factor (`--max-ice-factor 1.2`), run a short segment, `check`.

**Stage 2 — production, one case in its watch window at a time:**

```bash
# on the HPC, per case, at a segment boundary
./run_trend_batch.sh --cam --cice ... <case>             # exocam-trend: fresh series
exocam-accelerate advise <trend_dir> <case> --json advice.json
exocam-accelerate taper --advice advice.json --archive <DOUT_S_ROOT> \
    --json tapered.json                                   # + tapered.taper.nc

# stop the chain first: no job for the case may be queued or running
exocam-accelerate jump --rundir <rundir> --advice tapered.json \
    --archive <DOUT_S_ROOT> --dry-run                     # pre-flight + preview
exocam-accelerate jump --rundir <rundir> --advice tapered.json --archive <DOUT_S_ROOT>

# resubmit a SHORT continuation segment (e.g. STOP_N=5 years), regenerate
# trends when it ends, then:
exocam-accelerate check <trend_dir> <case> --rundir <rundir>
#   exit 0  PASS  -> resume normal segments
#   exit 10 WAIT  -> another short segment, check again
#   exit 20 FAIL  -> roll back and resubmit:
exocam-accelerate rollback --rundir <rundir> --archive <DOUT_S_ROOT>
```

`jump` refuses (exit 3) unless: no SLURM job named after the case is queued or
running; every `rpointer.*` points at the restart being edited; no history in
the run directory is past that date; a verified rollback source exists — the
archived restart set for that date, pristine, via `--archive` (mandatory; the
only escape hatch is `--allow-no-archive-rollback`, which instead verifies
every `rpointer.*`-named component restart is still present in the run
directory, and is only a valid rollback source until the next `st_archive`
sweep); and the advice is bound to this jump — a known `schema_version`, the
same case, and at most 5 years old. It keeps the pristine `cice.r` and the
jump log in `run/exocam_accelerate/` (not next to the restart — CESM's
`st_archive` sweeps `${CASE}.cice.r.*`).

`advise` refuses (no ice factor; reasons explain why) rather than merely
warning when a diagnostic needed to trust the jump is missing or fails its
gate: no `ICEFRAC` series (settled-edge check), no `qi` series or `qi`/`hi`
drifting more than `--max-enthalpy-drift`, a rejected `TS(N)`/`Tsfc(N)`
reference, or an explicit `--n-target` that does not lie strictly between the
current imbalance and the fitted asymptote. Saved advice JSON carries a
`schema_version` and a `provenance` map (input trend file name → sha256),
binding it to the exact series it was fit from.

`check` scores the post-jump run against the relations the advice was fitted
on, evaluated at the ice the run actually has: the jump landed (first
post-jump `hi` ≈ expected), `N` sits on the conduction law, `TS` on its linear
`TS(N)` relation, and the ice is not melting back — after 2 adjustment years,
PASS needs 3 settled years. If the jump log's `TS(N)` reference is missing or
was not accepted by the advisor, `check` cannot PASS under the standard
`aqua_ice` protocol — once enough settled years have accumulated it returns
FAIL (not WAIT: no amount of additional settled data fixes a reference that
was never accepted). On null jumps in the 15 grp3 runs it false-alarms in
~1 % outside pt10 (whose own regime shifts it flags).

`rollback` restores the archived restart set for the jump date (all
components and rpointers — the run directory copy was the one edited), saves
the current rpointers, and lists the discarded segment's output without
deleting it: move that aside before any trend analysis.

The next `advise` detects the jump (a one-year step in `hi`) and fits
post-jump native annual means only; it needs 2 + 15 settled years after a jump,
so one cycle is ~20 model years for ~180 skipped. Defaults: remove half the
current imbalance (`--n-fraction 0.5`), target clipped to 5x the window's
N-range, ice factor capped at 1.5 (`--max-ice-factor`); hard bound 2.0 whatever
the flags say. Snow is left alone unless `--snow-factor` is given.

## Jumping a hot, ice-free case (`advise-ocean`)

For slab-ocean runs with no sea ice whose slow drift is the surface
temperature (e.g. thick CO2 atmospheres at 340–375 K). The jump shifts the
slab temperature `somtp` in `docn.r` (the file `rpointer.ocn` names), sized
from a Gregory line `TS = c0 + c1·N` fitted against the **surface** imbalance
`energy_bot`: that is the slab's own equilibrium condition, whereas
`energy_top` also carries the heat the atmosphere stores while it warms —
mostly the latent heat of the vapor column, 2–3× the ocean's heat capacity at
365–375 K. The somtp increment is scaled by the measured heat ratio
C_total/C_ocean so that TS lands on target after the atmosphere takes its
share. Details and the evidence: `docs/ocean-jump.md`.

```
exocam-accelerate advise-ocean TRENDDIR CASE --json adv.json   # --rundir RUN picks up a past jump
exocam-accelerate pattern --advice adv.json --archive ARCH \
    --domain-file DOMAIN --json adv_p.json                     # optional: measured warming pattern
exocam-accelerate jump --rundir RUN --advice adv.json --archive ARCH --dry-run
exocam-accelerate jump --rundir RUN --advice adv.json --archive ARCH
# resubmit a short continuation, regenerate trends, then:
exocam-accelerate check TRENDDIR CASE --rundir RUN
exocam-accelerate somtp-map --rundir RUN -o somtp.nc          # somtp on the lat-lon grid
```

`--probe` (or `--probe-dt K`) instead steps TS ahead by the recent trend ×
`--probe-years` to *map* the trajectory where the Gregory line is not yet
constrained; `check` then measures λ and the implied equilibrium from the
probe's lever arm, and fails a probe after which the imbalance grows.

**Atmosphere too (`atm-profile`).** At 365–375 K an ocean-only jump is mostly
drained into the vapor column. `atm-profile --advice ADV --archive ARCH --json
C.json` couples the advice with a cam.r edit: T per level × a profile measured
from two archived cam.i, q at fixed RH, dry mass per layer kept (PS rises),
`TEOUT` recomputed so CAM's energy fixer does not undo it, and the stratiform
scheme's previous-step T/q shifted. `jump --advice C.json` then edits docn.r
and cam.r together (model constants from the run's `atm.log`); `check` and
`rollback` treat them as one jump.

The Gregory advisor refuses when the run still has ice, when TS does not rise as the
imbalance falls, when the current state is off the fitted line (a steepening
relation — sensitivity rising with temperature), or when the jump would be
smaller than the year-to-year TS scatter. The domain file is read from the run
directory's `docn_ocn_in` when `--domain-file` is not given. `--heat-ratio`
overrides the measured ratio; `check` reports the ratio the first post-jump
year implies.

## Viewing a case (`exocam-accelerate view`)

```
exocam-accelerate view DIR          # -> http://127.0.0.1:8765 (next free port if busy)
```

`DIR` holds exocam-trend `<case>_*_{cam,cice}.txt` files and, optionally, jump
logs (`*.nc.accel.json`, copied from `run/exocam_accelerate/`; one level of
subdirectories is searched). The page shows, per case: the advice (or refusal
reasons), the phase-space picture N vs hi with the fitted conduction law and
the jump's landing point, TS vs N, the model years a jump skips (Stefan
growth on the model clock, post-jump data in equivalent years), time series
with fit window and predicted levels, and each jump's `check` verdict.
Phase-space panels are drawn as Gregory plots (imbalance on y, temperature
on x), and every time-series row can show its time derivative (values / d/dt /
both; centered 5-yr slope, piecewise between jumps). Hot, ice-free cases get
the ocean layout instead: `energy_bot` against TS with the Gregory line (and every fit window tried), the same years against
`energy_top` (the offset is atmospheric storage), `energy_bot` vs `energy_top`, the
one-box projection a jump short-cuts, and lat-lon maps of any `*.latlon.nc`
(`somtp-map`) or `*.pattern.nc` (`pattern`) files in `DIR`.
Controls re-run the advisor live (what-if); nothing is written. The page
refreshes when files in `DIR` change. Fetching new trend data from the HPC is
a separate step.

## Related tools

- [`exocam-trend`](../exocam-trend) (dependency) — produces global-mean
  time-series diagnostics from monthly model output; its text output is one
  input pathway for tendency measurement (`trend_io`).
- [`exocam-casemgr`](../exocam-casemgr) (reference only) — case build/run/data
  management. This tool will eventually be invoked from casemgr as an
  external, but does not import from it.

## Install (development)

```bash
pip install -e ".[dev,netcdf]"   # netcdf needed only for `jump`/`restore`
pytest
```

## Status

Experimental. Pure-computation safeguards, the Tier-0 offline hindcast
harness, the phase-space jump advisor, and the `aqua_ice` restart writer
(cice.r, in place with backup, continuation), the post-jump `check`, pre-flight
safety checks and whole-set `rollback` are implemented and tested (152 unit
tests at the time). The advisor is validated offline on 15 cold cases, and the
first real ice jumps (grp4 pt01/pt03, 2026-09-28) passed. The `som_ocean`
(docn.r somtp) jump for hot, ice-free runs is implemented and validated
offline on atlasfu D1–D5; no real ocean jump has been made yet. Run
orchestration (submitting and polling segments), land and atmosphere plugins
are not implemented.

**Use with care and caution — model behavior is not always straightforward.**

## Glossary

Plain-language definitions of the terms used above and in the CLI output. Sections of this README and `docs/` carry the detail; entries here only say what a word means.

### Workflow and files

- **jump**: one edit of a restart file that moves the model state a number of years ahead along its own trajectory toward equilibrium (`jump`). Followed by a short run to see whether the model agrees (`check`).
- **advice / advisor**: `advise` (ice), `advise-ocean` (hot, ice-free) and `taper` read exocam-trend series and propose a jump (ice factor or somtp increment, the target, the model years skipped) or refuse with reasons. `--json` saves the advice, which `jump --advice` consumes; it carries a `schema_version` and the sha256 of its input series.
- **continuation**: the restarted run is a plain continuation (no hybrid or branch reboot). The edited restart is read as if the model had written it, so nothing in the case setup changes.
- **in place**: the restart file in the run directory is edited itself, after a pristine copy has been saved. No new file name enters `rpointer.*`.
- **pristine backup (`.pre-accel.nc`) and jump log (`.accel.json`)**: copy of the unedited restart and the record of what was changed (factor or increment, restart date, the advice). Both live in `run/exocam_accelerate/`, not beside the restart, because `st_archive` sweeps `${CASE}.cice.r.*`.
- **pre-flight**: the checks `jump` runs before editing (`--dry-run` stops after them): no queued or running job, `rpointer.*` consistent, no history past the restart, a verified rollback source, advice matching the case and recent.
- **rollback / restore**: `rollback` resets the whole restart set (all components and rpointers) to the archived pre-jump set. `restore` puts back just the pristine file (`--ocean` for docn.r) before resubmitting.
- **`cice.r`**: CICE sea-ice restart: ice volume (`vicen`), enthalpy (`eicen`), area (`aicen`) per category.
- **`docn.r`**: data-ocean (slab) restart; holds `somtp`, the mixed-layer temperature in K on the docn grid (`gsize = ni*nj`, longitude fastest).
- **`cam.r`**: CAM atmosphere restart: T, q, winds, `DELP`, `PS`, `TEOUT` and the previous-step fields.
- **`cam.i`**: CAM initial-condition file written at the end of a year in the archive. Two of them, years apart, give the measured warming profile (`atm-profile`).

### Energy balance

- **`energy_top`** (N_TOA, TOA balance): net top-of-atmosphere flux. The convergence criterion, and the default imbalance for `advise-ocean` (`--imbalance`). N = 0 means equilibrium.
- **`energy_bot`**: net flux into the slab ocean. About 3 times noisier year to year; kept for information and selectable with `--imbalance energy_bot`.
- **TOA-surface gap** (`energy_top - energy_bot`): heat the atmosphere stores while it warms, mostly the latent heat of the water-vapor column. It is not a leak and closes at equilibrium (`docs/ocean-jump.md`).
- **N now** (`N_now`, `N_now_se`, `N_last_year`): the imbalance now, as the mean of the last 5 native annual values with its standard error. `N_last_year` is the single latest year, which can swing by several W/m2.
- **`N_target`, `n_fraction`**: the imbalance the jump aims for, by default `N_now*(1 - n_fraction)` with `--n-fraction 0.5`, or an explicit `--n-target`.

### Heat capacity

- **C_ocean, C_atm, C_total**: heat taken up per K of surface warming, in W yr m-2 K-1: the slab's, the atmosphere's (in hot cases mostly vapor latent heat, rising steeply with T), and their sum. Measured from N against TS over the fit window.
- **ocean heat fraction** (`f_ocean = C_ocean/C_total`; `--ocean-fraction`, JSON `heat.ocean_fraction`): the share of a TS change's heat the slab holds. The somtp step is the TS step divided by `f_ocean`. The **atmosphere heat fraction** is `1 - f_ocean` (`heat.atm_fraction`).
- **heat ratio** (`heat_ratio`, `--heat-ratio`): the old name, `1/f_ocean`. Still written to the JSON and still accepted as a hidden flag; `check` reports `heat_ratio_implied`.
- **effective heat capacity** (`C_eff = N / (dTS/dt)`): what the run actually shows over the window, also given as metres of sea water for comparison with `hblt`.
- **relaxation time** (`tau = C/alpha`): e-folding time of the imbalance in the one-box model, `C dTS/dt = N`, `N = alpha (TS_eq - TS)`.
- **years skipped**: model years the jump saves, `tau * ln(N_now/N_after)` for the ocean jump (Stefan growth for the ice jump).

### Phase space and feedback

- **Gregory analysis / phase space**: plot N against TS (or ice thickness). The state moves along a curve toward N = 0, and a jump moves it along that curve instead of waiting. The viewer draws N on y, T on x.
- **feedback parameter** (`lambda`, `alpha`): `-dN/dTS`, in W/m2/K. Positive means restoring.
- **differential feedback parameter** (`alpha_diff = -dN/dTS` at the current state; Gregory et al. 2004): the local slope of N(TS). For hot planets the feedback depends on T and changes regime (Wolf et al. 2018, JGR-A), so one line over the whole history misleads.
- **local curve** (`--fit local`, default): a quadratic N(TS) about the current state on native annual means. Curvature is used only when resolved at 2 sigma, otherwise the tangent. Refuses if `alpha_diff <= 0` or not resolved at 2 sigma, or if the curve turns over before N reaches the target.
- **line** (`--fit line`): the original single Gregory line `TS = c0 + c1*N` per window. Checked by "on the line": the current state must lie within `max_offline` of it.
- **equilibrium**: N = 0; TS there is `TS_eq`.

### Trend data

- **exocam-trend `data/*.txt`**: whitespace-separated monthly series per variable, read by `trend_io`. Each variable has `native`, `int1` and `int2` columns: the monthly values, a short running mean, a long running mean. `--which` picks one (default `int2`); the local fit uses native annual means so years are independent.

### Safeguards

- **trustworthiness gate** (time-domain): refuses a step when the trend window shows accelerating curvature (the drift is speeding up, near a threshold) or a significant sign reversal of the tendency (2 sigma). Decelerating curvature, the normal shape of convergence, passes.
- **hard clip** (`--max-dt`, `--max-ice-factor`): absolute cap on the step, whatever the fit says.
- **noise floor**: a jump smaller than the interannual TS scatter is refused as not worth an edit.
- **extrapolation limit** (`--max-extrapolation-ratio`, 5): the target may not lie more than that multiple of the window's N-range beyond the data.
- **`--override-gate`**: turns the overridable refusals (feedback unresolved, state off the curve, curve turning over, low correlation, the time-domain gate for a probe) into warnings recorded as `gate_overridden`. `alpha_diff <= 0`, sea ice present and missing data still refuse.

### Probe and post-jump check

- **probe** (`--probe`, `--probe-dt K`, `--probe-years`): a deliberate perturbation to measure the response, not sized to equilibrium. The step is the recent TS trend times `--probe-years` (default 15), or an explicit K.
- **lever arm**: the TS difference between the probe's settled years and the years before it, which gives the response a baseline that the undisturbed run lacks.
- **post-jump check** (`check`): verdict **PASS** (exit 0), **WAIT** (10, run more) or **FAIL** (20, roll back). For a probe FAIL is mechanical only: it did not land, or the data are non-finite.
- **landed**: the first post-jump year moved by the amount the edit implies (`hi` for ice, TS for somtp). **Settle years** (default 2) are skipped after a jump, and PASS needs `min_years` settled years (3 by default; 5 for a probe).
- **energy readout**: lambda and the implied equilibrium from settled post-probe N against the years before. Uses `energy_top`; `energy_bot` is reported for information only.
- **TS-trajectory readout**: TS drift rate before versus after the probe gives `tau` and `TS_eq`. Much quieter than the energy readout on noisy runs.
- **above / bracketed**: whether the equilibrium lies beyond the probe level (hotter still) or between the pre- and post-probe states (the probe overshot). Readouts that disagree give WAIT.
- **runaway greenhouse suspected**: a warning when a warm probe is followed by no restoring response at 3 sigma (N grew, or the TS drift rose). Verdict WAIT with no rollback recommendation; the call is yours.
- **implied ocean heat fraction**: what the first post-jump year says `f_ocean` was; use it for the next jump's `--ocean-fraction`.

### Coupled atmosphere jump

- **coupled jump** (`atm-profile`): edits `docn.r` and `cam.r` together, with the atmosphere following the ocean. `check` and `rollback` treat them as one jump.
- **measured profile**: horizontal-mean warming per level between two archived `cam.i`, per K of surface warming. T is scaled by it.
- **q at fixed RH**: humidity is raised with T so relative humidity is unchanged; warming T alone leaves the air undersaturated and the ocean re-evaporates.
- **dry mass kept**: each layer's dry air mass is conserved, so adding vapor raises `PS`.
- **`TEOUT`**: column energy at the end of the last physics step. It is recomputed because CAM's energy fixer would otherwise remove the jump on the first step.
- **troposphere only**: nothing is changed at p < 100 hPa (`--ceiling`); the jump tapers linearly in log p from full strength at 200 hPa (`--taper-bottom`).

### Ice jump

- **ice factor**: one multiplier on `vicen` and `eicen`; `aicen` is untouched. Capped at 1.5 by `--max-ice-factor`, hard bound 2.
- **Stefan law** (`N = a + b/hi`): growth limited by conduction, so the imbalance falls as ice thickens. The target N sets the jump size.
- **tapered jump** (`taper`): the factor is applied per cell as `1 + (F-1)*weight`, with weight from the Stefan product `dh/dt*h` (ramp `--ramp 0.05 0.40`); partial-ice cells are not scaled.
- **effective factor** (`effective_factor`): area-mean ice after over before, what the tapered jump amounts to globally. `check` scores against the law with it.

### Pattern

- **pattern** (`pattern`): per-cell somtp weights from the measured warming pattern (somtp from two archived `docn.r`), so the increment follows where the planet warms rather than being uniform. Area mean stays the sized increment; `--max-weight` limits any one cell.
