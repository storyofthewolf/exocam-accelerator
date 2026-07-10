# Claude Code Prompt: `exocam-accelerate`

> **Before running:** fill in the two sibling repo paths, and attach in the
> session: the `accelerate.py` prototype, the two precedent PDFs (Wordsworth
> et al. 2013; Turbet et al. 2021), and the reference restart-file header dumps.

---

## Project: `exocam-accelerate` — a new standalone repository

You are starting a new Python repository, `exocam-accelerate`. This is a
scientific tool to accelerate the convergence of ExoCAM/CESM climate
simulations by extrapolating slow-evolving restart-file fields toward
equilibrium, based on measured time-series tendencies. It is a generalization
and safe redesign of an existing naive prototype script (`accelerate.py`),
which the user will provide.

**Do not implement the CESM restart-write / run-continuation layer in this
session.** This session is scaffolding plus investigation only. The algorithm
and its safeguards are designed; the CESM-integration mechanics are an open
question to be investigated against real code, not guessed at. See phasing
below.

---

## Scientific method (settled — implement the pure-computation parts of this)

The core operation is a Wordsworth-style forward-Euler extrapolation of a slow
state variable:

```
X_new = X + <dX/dt> * Δt
```

where `<dX/dt>` is a tendency measured from a simulation's time-series, and
`Δt` is a multi-year acceleration timestep. Precedent: Wordsworth et al. 2013
(Icarus, early Mars ice equilibration, §2.3) and Turbet et al. 2021
(day-night cloud asymmetry / early Venus, Methods §4). The user can supply
both.

The reusable core is the set of safeguards, not the per-variable scaling.
Implement these as pure functions operating on trend data and arrays (no file
I/O dependencies):

1. **Horizontally-averaged tendency.** Tendency should be taken from the
   horizontally-averaged per-layer trend, not pointwise gridcell tendencies
   (pointwise amplifies dynamical noise; the horizontal mean isolates the slow
   coherent drift toward equilibrium).

2. **Trustworthiness gate before any step.** Check that the second derivative
   (curvature, d²X/dt²) over the trend window is small and the tendency sign is
   stable. Large curvature signals proximity to a nonlinear climate-feedback
   threshold (ice-albedo, runaway greenhouse), where linear extrapolation is
   unsafe — the correct action there is to **refuse** the step, not shrink it
   blindly.

3. **Hard per-step magnitude clip** independent of what the extrapolation
   suggests (Turbet clips ΔT at 50 K per step for numerical stability).
   Configurable per variable.

4. **Decreasing Δt schedule.** Both precedents ramp Δt down as convergence
   approaches — large jumps early, small jumps near equilibrium.

5. **Post-jump physical-constraint re-imposition, per variable** (Wordsworth:
   enforce ice ≥ 0, conserve total water mass after each jump). This is the
   only variable-specific piece and should be pluggable.

6. **Post-step consistency check (design only this session).** After the model
   runs forward from an accelerated restart, verify the model's own tendency
   continued along the predicted trajectory rather than correcting away from
   it; if inconsistent, the step crossed a nonlinearity and must be rolled
   back. Scaffold the interface for this now; don't implement the run-forward
   half.

---

## Dependency boundary (important)

- **`exocam-trend`** (sibling repo, path: ../exocam-trend) is a **dependency**: this
  tool consumes a series of monthly mean climate output files and creates time-series outputs. Read its public interface and evaluate if this is resuable for import from it, or if duplicative methods should be created in this repo.  Note, exocam-trend currently only produces time-series of global mean quantities, no spatial differentiations.

- **`exocam-casemgr`** (sibling repo, path: ../exocam-casemgr) is **reference-only** for
  now: read its restart-file handling and run-continuation / hybrid-run setup
  logic to understand the mechanics, but do NOT import from it or couple to its
  internals. This tool will be wired into casemgr as an external later; keep
  the boundary clean.

Short version: **depend on exocam-trend, learn from exocam-casemgr, import from
neither casemgr's guts.**

---

## First-session tasks

1. **Scaffold the repo:** standard Python package layout, `pyproject.toml`, a
   `CLAUDE.md` capturing this design and the settled-vs-open distinction, tests
   dir, README stating the scientific purpose and citing the two precedents.

2. **Implement the pure-computation safeguard functions** (items 1–4 above)
   with unit tests, operating on arrays / trend structures — no file I/O, no
   CESM coupling.

3. **Design (interfaces only, no implementation)** the variable-plugin
   structure for item 5 and the consistency-check hook for item 6.

4. **Investigate and write up, do not implement:** read `exocam-casemgr`'s
   restart handling and the CESM/ExoCAM run-continuation code to characterize
   the central open question — when an accelerated restart should be applied by
   modifying `.r.` files for a simple continuation run versus modifying `.i.`
   files and rebooting as a hybrid run, and how each path treats the perturbed
   fields (which fields are read verbatim on continuation vs. re-derived on
   hybrid init). Produce a findings document in the repo
   (`docs/restart-integration-questions.md`) laying out what you found and what
   still needs the user's decision. This is the design input for the next
   session.

Read the prototype `accelerate.py`, the two precedent papers, and the
`exocam-trend` interface before writing code. Ask the user for the sibling repo
paths and any files not provided.

---

## Reference restart files (metadata-only at this stage)

The user will provide a set of ExoCAM restart files from a real simulation
checkpoint, to anchor the task-4 investigation against ground truth rather than
only inferring file layouts from code. **Inspect these with header dumps only
(`ncdump -h`, or equivalent netCDF metadata reads). Do NOT read the full data
payloads — these files are hundreds of MB and only their structure is needed
now.** Reading array values is explicitly out of scope for this session.

Restart file examples are found in
- /restart_file_examples/cam_aqua_fv: aquaplanet (ocean-only) model configuration restart file set
- /restart_file_examples/cam_land_fv: landplanet (land-only) model configuration restart file set
- /restart_file_examples/cam_aqua_fv: mixed (ocean + land)  model configuration restart file set

From the headers, characterize for the findings document:

- Which variables and dimensions live in each file type (`cam.r`, `cam.rs`,
  `cice.r`, `clm2.r`, `docn.r`, and the `.i.` initial file if provided), so the
  acceleration targets can be mapped to their actual file homes.

- The structural difference between `.r.` (restart) and `.i.` (initial) files
  where both are provided — which fields appear in each. This directly informs
  the `.r.`-continuation vs `.i.`-hybrid-reboot decision that is the session's
  central open question.

Note in the findings document if the `.i.` files were not available in the
reference set (this itself is informative — it indicates the automated workflow
may need to regenerate them). The reference set should ideally come from a
cold/frozen-regime run, since the aqua_ice snowball case is the first
acceleration target and its ice-energy fields (`eicen`, `vicen`, `esnon`,
`vsnon` in `cice.r`) are the variables the first plugin will touch; note the
regime of the provided set regardless.

---

## Potential Restart Methods
1. continuation
reads rpointer.* files, reads the files specified in the rpointer.* files.
intent is to restart the model bit-for-bit from a restart file set.
2. hybrid
mixes restart files (.r.) with initial conditions files (.i)
cam.i. file is used
ocean initial file is also used (pop_frc.*) found in /ocean_initial_examples
pop_frc.4x5d.090130_aquaplanet_300Kiso.nc is used for aquaplanets
pop_frc.gx3v7.110128_annual_mean.nc is used for mixed cases with Earth continents.
3. new "initial" simulation
- start a new simulation using only initial condition files produced by a previous simulation



## Note for later development (do not implement)

A second acceleration mechanism, Turbet-style in-situ radiative heating-rate
multiplication (multiply deep-layer heating rates by a factor, ramp
0.5 → 0.3 → 0 as convergence approaches), is deferred. It requires modifying
the GCM Fortran source rather than manipulating restart files, so it is out of
scope for this file-manipulation-based tool for now. Recorded here so the
design space is not forgotten.


