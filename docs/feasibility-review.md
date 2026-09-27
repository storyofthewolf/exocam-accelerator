# Feasibility review and recommended validation plan

Review date: 2026-09-27

Reviewed state: `feature/advise-and-jump` at `b75fa23`

## Decision

The core idea is scientifically plausible for a narrow target: cold
slab-ocean aquaplanets with a settled ice edge whose slow equilibration mode is
conduction-limited sea-ice thickening. The phase-space formulation
`N = a + b/hi` is a much stronger basis for this target than linear
extrapolation in time.

The project is ready for a controlled canary after the fail-open safety issues
below are corrected. It is not yet validated for unattended production use or
as a general ExoCAM convergence accelerator.

The current evidence supports an experimental `aqua_ice` state-jump tool. It
does not yet support claims about land ice, mixed-surface cases, moving ice
edges, deep-atmosphere equilibration, or climates near an ice-albedo
bifurcation.

## What is already strong

- The project abandoned the weak original assumption that the slow state can
  be advanced reliably with a linear tangent in time. The Tier-0 results show
  systematic overshoot for sea-ice thickness and motivate the phase-space
  formulation.
- The empirical Stefan relation has a physical interpretation and a strong
  within-suite fit. The recorded analysis reports `r^2 >= 0.98` in 12 of 15
  cold cases and `0.91--0.98` in the remaining cases.
- Scaling `vicen` and `eicen` together preserves enthalpy per unit ice volume.
  Leaving `aicen` unchanged also avoids an immediate imposed change in ice
  fraction and albedo.
- Restart edits are recoverable in the intended archived workflow, jumps are
  reapplied from a pristine copy rather than compounded, and the code records
  advice and provenance.
- The code is compact and well separated into pure fitting logic, restart I/O,
  run-state checks, and CLI behavior. At the reviewed state, all 152 unit tests
  pass.

## Findings that should be resolved before a live jump

### 1. Require a verified full restart set for rollback

The README says that `jump` refuses unless a pristine archived restart set
exists, but `preflight()` currently treats a missing `--archive` as a warning.
After the next archive cycle, the run directory may no longer contain the
other components at the pre-jump date, making whole-set rollback impossible.

Production behavior should fail closed unless either:

1. the complete archived restart set has been located and verified; or
2. a separate no-archive mode has verified that every component restart and
   rpointer needed for rollback is retained locally.

### 2. Make the scientific gates fail closed

The implemented jump assumes a settled ice edge, yet missing `ICEFRAC` is only
a warning and still allows advice with an ice factor. A rejected `TS(N)` fit,
including disagreement that may indicate a feedback kink, is also only a
warning. These conditions should prevent a jump for the initial production
workflow.

The post-jump checker should also honor the recorded `accepted` flag for the
temperature relation. It currently uses a stored linear fit even when the
advisor rejected that relation. If the required temperature reference is
missing or rejected, the checker should not return `PASS` under the standard
`aqua_ice` protocol.

The minimum fail-closed inputs should be:

- accepted hyperbolic `N(hi)` fit;
- present and settled `ICEFRAC` series;
- present `qi` series with acceptably stable `qi/hi`;
- accepted `TS(N)` reference;
- current advice for the exact case and restart date; and
- a verified rollback source.

### 3. Validate the spatial and category state

The advisor derives one factor from global-mean `hi` and applies it uniformly
to every grid cell and thickness category. Natural conductive growth is
approximately additive in `h^2`, so thin and thick ice do not generally evolve
by the same multiplicative factor. A stable global ice fraction does not prove
that the regional thickness pattern or ice-thickness distribution remains
valid after uniform scaling.

CICE's first-step remapping may reconcile category bounds, but that response is
part of the scientific experiment and cannot be assumed. A global-mean
post-jump check can pass while regional ice thermodynamics or category
populations are distorted.

The canary must therefore inspect maps and category-resolved fields in addition
to global `N`, `TS`, and `hi`.

### 4. Validate the intervention, not only its destination

The offline advisor hindcast shows that a recommended state lies near a state
the uninterrupted model later reaches. It does not show that writing that
state into `cice.r`, while leaving the atmosphere, slab ocean, snow, and
coupler restart state unchanged, stays on the same attractor.

The reported 165 origins are repeated windows from 15 related cold-start
trajectories. They are valuable calibration samples but not independent
climate realizations. Generalization should be evaluated by holding out entire
runs and, later, distinct regimes.

### 5. Bring the scientific evidence into the repository

The phase-space bake-off, advisor hindcast, and false-alarm drivers are
currently recorded as scripts outside the repository. The source time-series,
case manifest, exact smoothing configuration, and machine-readable result
tables are also absent. The unit tests verify the algorithms on synthetic
trajectories but cannot reproduce the headline scientific numbers.

Record enough of the workflow to regenerate the reported tables and figures:

- analysis drivers;
- a case manifest with case descriptors and data provenance;
- exact `exocam-trend` command and smoothing windows;
- machine-readable summary tables; and
- either a small distributable fixture or a documented external-data manifest
  with checksums.

### 6. Validate all numerical controls

The CLI should reject invalid domains before producing advice. At minimum:

- `0 < n_fraction < 1` for the standard workflow;
- `1 <= max_ice_factor <= 2` for a thickening-only cold-case jump;
- `max_extrapolation_ratio > 0`;
- positive window and post-jump durations; and
- an explicit `N_target` that advances toward the permitted side of the fitted
  asymptote without reversing the intended equilibration direction.

Advanced experiments that intentionally thin ice or cross these bounds should
use a separately named override rather than overloading the normal path.

### 7. Reconcile the documentation and active test plan

The README still says that the constraint plugin, post-step check, and restart
writer are unimplemented before later describing them as implemented. The
older Tier-2 plan also specifies a combined `cice.r + somtp` perturbation and a
continuation/hybrid comparison, while the implemented decision is a cice-only
continuation.

The active test plan should describe the exact deployed intervention. Hybrid
should remain a fallback experiment if the simpler continuation exhibits a
persistent coupling shock.

## Recommended validation sequence

### Stage 0: close the pre-flight gaps

Before a live restart is edited:

1. make archive/full-set rollback verification mandatory;
2. make missing or rejected scientific diagnostics block advice;
3. validate CLI parameter domains;
4. make `check` honor fit acceptance;
5. add an advice schema version and bind advice to the case, restart date, and
   relevant input-series provenance; and
6. make the offline scientific results reproducible from a recorded manifest.

### Stage 1: sham restart test

Create two identical clones from one restart set:

- control: untouched continuation;
- sham: pass `cice.r` through the writer with factor 1.

Compare model field values after the first timestep, first day, first month,
and short segment. The test should distinguish harmless NetCDF metadata changes
from model-state differences. Verify all component restarts, rpointers, logs,
and rollback behavior.

Acceptance criterion: the sham and control are bit-for-bit identical in model
fields, or any difference is understood, bounded at roundoff, and does not grow.

### Stage 2: gentle cice-only canary

Use one clone pair and apply a factor of approximately `1.05--1.10` to the
canary. Do not add `somtp` or switch to hybrid in this first intervention.

Check at high cadence during the initial adjustment, then over several years:

- global `N`, `TS`, `ICEFRAC`, `hi`, and surface energy fluxes;
- maps of `hi`, `Tsfc`, ice fraction, and conductive/basal fluxes;
- per-category `aicen`, `vicen`, and `eicen`;
- thickness-distribution and enthalpy-per-volume changes;
- CICE remapping during the first steps;
- atmosphere/ice/coupler log warnings; and
- whether rollback fully restores the control starting state.

Acceptance criterion: no numerical shock, no persistent regional or category
artifact, and adjustment toward the predicted `N(hi)` and `TS(N)` relations.

If continuation exhibits a persistent coupler artifact, repeat the same
perturbation as a controlled continuation-versus-hybrid comparison. Do not
introduce hybrid complexity unless this evidence requires it.

### Stage 3: causal twin-convergence experiment

Restart accelerated and untouched twins from an early checkpoint in several
held-out cases. Include at least:

- a cold waterbelt with a settled edge;
- a full snowball;
- a case near the edge of the advisor's acceptance region; and
- a case not used to tune the default thresholds.

Declare the equilibrium and equivalence criteria before running. Compare:

- final global energy balance and drift;
- global and spatial temperature;
- sea-ice area, volume, snow, and category distributions;
- relevant surface and TOA budgets;
- internal variability over a common averaging interval; and
- simulated years, wall time, and core-hours to convergence.

Acceptance criterion: the accelerated twins reach the same equilibrium basin
and statistically equivalent climate while delivering a material reduction in
core-hours. This is the experiment that can support the claimed acceleration,
including the estimate that one capped jump may replace roughly 180 simulated
years in the existing cold-case suite.

### Stage 4: production integration

Only after the twin experiment succeeds should the workflow be connected to
`exocam-casemgr` for submission, monitoring, resumable state, and archival.
Keep one case in a watch window at a time until several independent jumps have
passed.

Defer learned parameters, `somtp`, land, and atmospheric plugins until the
`aqua_ice` intervention has demonstrated end-to-end benefit. Each new problem
class needs its own physical fit, constraints, diagnostics, and validation
suite.

## Go/no-go gates

| Gate | Required result |
|---|---|
| Code and pre-flight | Fail closed on missing rollback or diagnostics |
| Sham restart | No growing difference from untouched continuation |
| Gentle canary | Stable coupled response and valid spatial/category state |
| Twin convergence | Same equilibrium basin within internal variability |
| Performance | Material measured core-hour reduction |
| Generalization | Success on held-out runs and more than one cold regime |

The immediate recommendation is **go for Stage 0 and then one controlled
canary; no-go for production deployment or broader ExoCAM claims until the
twin-convergence gate passes**.
