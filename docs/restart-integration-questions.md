# Restart integration: findings and open questions

Session-1 investigation (2026-07-10) for the central open design question:

> When should an accelerated state be applied by **modifying `.r.` restart
> files and continuing** the run, versus **modifying/creating `.i.` initial
> files and rebooting as a hybrid run** — and how does each path treat the
> perturbed fields (read verbatim vs. re-derived at init)?

Sources: `ncdump -h` header dumps of the three reference restart sets under
`restart_set_examples/` (metadata only; no data payloads read),
`exocam-casemgr`'s branch/hybrid build logic and `runmgr.py` continuation
logic (reference only), the `rpointer.*` files, and the two precedent papers.

**Nothing here is implemented.** This document is the design input for the
next session's decisions.

---

## 1. The three restart paths (mechanics as found)

### 1a. Continuation (`CONTINUE_RUN=TRUE`)

`runmgr.py continue` sets `CONTINUE_RUN=TRUE` and resubmits. The driver reads
`rpointer.drv`; each component reads its own `rpointer.*` naming its restart
file. `rpointer.atm` names the master `cam.r.` file and lists the companion
files (`cpl.r`, `cice.r`, `docn.r`, `docn.rs1.bin`, or `clm2.r`). Intent is a
**bit-for-bit** resume: every field in every `.r.` file is read verbatim;
nothing is re-derived. An accelerated state applied this way must therefore be
internally consistent across *all* the coupled state the files carry (see §4).

### 1b. Hybrid reboot (`RUN_TYPE=hybrid`)

casemgr's `_build_branch_post_setup` shows the recipe: copy the reference
restart set into the new rundir, then
`xmlchange RUN_TYPE=hybrid, CONTINUE_RUN=FALSE, RUN_REFCASE=…, RUN_REFDATE=…`.
In a hybrid run:

- **CAM starts from the `.i.` file** (`ncdata`), *not* from `cam.r`. All the
  auxiliary physics-buffer state in `cam.r`/`cam.rs` is re-initialized.
- **The data ocean uses a `pop_frc.*` SOM forcing file**
  (`user_docn.streams.txt.som` in casemgr):
  `pop_frc.4x5d.090130_aquaplanet_300Kiso.nc` for aquaplanets,
  `pop_frc.gx3v7.110128_annual_mean.nc` for mixed/Earth-continent cases.
  Confirmed in `ocean_initial_examples/`: monthly climatology of `T`, `S`,
  `U`, `V`, `hblt`, `qdp` on the model grid.
- **CICE and CLM have no `.i.` files in the reference sets.** In CESM1-era
  hybrid runs these components still initialize from their `.r.` files from
  the reference set (CLM alternatively via `finidat`). So ice/land
  acceleration edits go into `.r.` files *on either path* — the hybrid/
  continuation choice only changes how the **atmosphere and coupler** state
  is treated.
- The model gets a fresh time axis and re-derives the coupler exchange state
  on the first coupling step.

### 1c. New "initial" run

Start a fresh case using only initial-condition files produced by a previous
simulation (`ncdata` = a `cam.i.`, plus `pop_frc`/`finidat` as applicable).
Cleanest state, most re-derivation, loses the most spun-up detail (no ice
state at all unless CICE gets a restart anyway — effectively the hybrid path
with extra steps for our purposes).

---

## 2. File inventory from the reference headers

All three sets are ExoCAM/CESM1.2.1, CAM4 FV 4°×5°, 72×46 grid; the
`cam_*_fv` sets use L51. `gsize = column = 3312 = 72×46`.

| File | Contents (acceleration-relevant) | aqua | land | mixed |
|---|---|---|---|---|
| `cam.r` | FV dycore prognostics `U V DELP PT Q CLDLIQ CLDICE PS` (lev,lat,lon); physics buffer (`QCWAT/LCWAT/TCWAT`, `CLD`, `QRS/QRL`, `tke/kvh/kvm`, `DTCORE`, `T_TTEND`, …); surface-flux snapshots; history-tape bookkeeping | ✓ | ✓ | ✓ |
| `cam.rs` | Coupler-exchange snapshot for the atm: `x2a_*` (incl. `Sf_ifrac`, `So_t`, `Si_snowh`, albedos) and `a2x_*` vectors | ✓ | ✓ | ✓ |
| `cam.i` | Clean prognostics on (time,lev,lat,lon): `T Q CLDLIQ CLDICE PS US VS` + surface `TS1–TS4 TSICE SICTHK SNOWHICE ICEFRAC` | ✓ | ✓ | ✓ |
| `cice.r` | Per-category ice state (ncat=5): `aicen vicen vsnon Tsfcn`; per-layer energies `eicen(ntilyr=20) esnon(ntslyr=5)`; dynamics/stress; `iceumask` | ✓ | — | ✓ |
| `clm2.r` | Column/pft-indexed land state: `T_SOISNO(column,levtot) T_GRND T_LAKE T_VEG T_REF2M*`; hydrology `H2OSOI_LIQ/ICE H2OSNO SNOWDP WA ZWT` | — | ✓ | ✓ |
| `docn.r` | Single vector `somtp(gsize)` — SOM prognostic temperature | ✓ | — | ✓ |
| `docn.rs1.bin` | Binary stream restart for docn forcing | ✓ | — | ✓ |
| `cpl.r` | Driver state: fraction fields, budgets, full `a2x/x2a/i2x/o2x/l2x` exchange vectors (154 variables in the aqua set) | ✓ | ✓ | ✓ |

**`.i.` availability:** `cam.i` **is present in all three sets** — the
automated workflow does not need to regenerate atmosphere initial files for
these cases (though a general tool must handle cases archived without
`inithist`). No `cice.i`, `clm2.i`, or `docn.i` exist; that is normal for
this model generation — hybrid runs swap only the CAM initial file and the
ocean forcing.

**Regime of the reference sets (matters for the first plugin):** none of the
three is a cold/frozen-regime run. The aqua set is THAI **Hab1**
(temperate aquaplanet, year 0217), the land set THAI **Ben1** (dry
landplanet, year 0259), the mixed set a modern-Earth-like case (year 0061).
Since `aqua_ice` (snowball) is the first acceleration target, a
frozen-regime restart set should be pulled before implementing the cice
plugin, so realistic `eicen`/`vicen` magnitudes and `iceumask` extent can
anchor the constraint logic.

---

## 3. `.r.` vs `.i.` — the structural difference (CAM)

The comparison is only possible for CAM (the only component with both):

- **`cam.i` carries temperature as plain `T`; `cam.r` does not have `T` at
  all.** The FV restart stores the dycore state — `PT` (potential
  temperature) and `DELP` (layer pressure thickness) — plus several
  *redundant copies* of thermal state in the physics buffer (`TCWAT`,
  `T_TTEND`, `DTCORE`, radiative heating rates `QRS/QRL`, …) and in the
  coupler snapshot (`cam.rs: x2a_Sx_t`, `a2x_Sa_tbot`, …). This is exactly
  why the prototype's `--atm_temp` mode found "no T or TS in cam.r".
- **Winds differ in representation**: `cam.i` has staggered `US(lev,slat,lon)`
  / `VS(lev,lat,slon)`; `cam.r` has A-grid-ish `U/V(lev,lat,lon)` plus all
  derived turbulence state (`tke`, `kvh`, `kvm`).
- **`cam.i` has no physics buffer and no flux snapshots** — on hybrid init
  the model re-derives all of it from the prognostics in the first steps.
- `cam.i` also carries simple surface fields (`TS1–4`, `TSICE`, `SICTHK`,
  `SNOWHICE`, `ICEFRAC`). With an active CICE + SOM configuration, at least
  the ice-related ones presumably come from `cice.r`/`docn.r` instead and may
  be ignored — **needs verification against CESM source** (open question O4).

**Implication:** perturbing atmospheric temperature through `cam.r` would
require consistently editing `PT` (with its `DELP` weighting) *and* every
redundant thermal copy in the physics buffer and `cam.rs`/`cpl.r` exchange
state — fragile and effectively a reimplementation of the model's own init.
Perturbing `T` in `cam.i` and rebooting hybrid lets CESM's own
initialization re-derive all secondary state. For the atmosphere, the hybrid
path is strongly indicated.

---

## 4. Per-target mapping: where each acceleration variable lives

| Target (prototype mode) | Variables | File | Path implication |
|---|---|---|---|
| Sea-ice energy/volume (`aqua_ice`) | `eicen vicen esnon vsnon` (+`aicen Tsfcn iceumask` for consistency) | `cice.r` only — no `.i.` exists | Must edit `cice.r` regardless of path. Open: does a *continuation* tolerate a perturbed `cice.r` when `cam.rs`/`cpl.r` still carry the old `Sf_ifrac`/`Si_snowh` exchange state, or does that require the hybrid reboot so the coupler re-derives fractions? |
| Atmosphere T (`atm_temp`) | `T` | `cam.i` (as plain T); only `PT`+pbuf copies in `cam.r` | Hybrid strongly indicated (§3). Per-layer horizontal-mean deltas map naturally onto `T(time,lev,lat,lon)`. |
| Land temperature (`land_temp`) | `T_SOISNO T_GRND T_LAKE T_VEG T_REF2M*` | `clm2.r` (column/pft vectors, `levtot=20` incl. 5 snow levels) | Edit `clm2.r` on either path. 1D column/pft indexing means the plugin needs the column→gridcell mapping for horizontal-mean tendencies; "levels" here are soil/snow depths, not atmosphere. |
| Land hydrology (`land_hydro`) | `H2OSOI_LIQ H2OSOI_ICE H2OSNO SNOWDP WA ZWT` | `clm2.r` | Same as above; water-mass conservation constraint spans several variables. |
| SOM ocean temperature | `somtp(gsize)` | `docn.r` | Editable in place for continuation; on hybrid the SOM is re-initialized from `pop_frc.*` — so ocean-temperature acceleration is *lost* on hybrid reboot unless a custom `pop_frc` is generated. Direct tension with the atmosphere's preference for hybrid. |

The coupled-consistency problem in one sentence: **continuation preserves
everything (including what we failed to update); hybrid re-derives everything
CAM-side (including what we might have wanted to keep, and it resets the SOM
ocean from a static forcing file).**

---

## 5. Open questions needing user decision / further verification

- **O1 — Path per target, or one path for all?** A mixed strategy looks
  natural from §4 (hybrid for atm-T; `.r.`-editing for ice/land/ocean), but a
  single acceleration step usually touches several targets at once (snowball:
  ice energy + surface/atm T + somtp). Do we standardize on one mechanism per
  *step*, and if so which wins when targets conflict (docn vs cam)?
- **O2 — Continuation after `cice.r` perturbation.** Does the CESM1 driver
  tolerate stale ice-fraction/exchange state in `cam.rs`/`cpl.r` on a
  continuation restart (re-converging within a coupling step or two), or does
  it produce shocks/crashes? Wordsworth and the prototype both effectively
  did in-place `.r.` edits and ran on — the prototype's empirical success at
  `aqua_ice` suggests continuation works for ice, but this session made no
  runs and the answer needs a controlled test on the target machine.
- **O3 — Who regenerates `.i.` files?** The reference sets have `cam.i`
  because these cases wrote `inithist`. For cases without one, the hybrid
  path needs a `cam.i` manufactured either by (a) enabling `inithist` and
  running one more segment, or (b) writing one from scratch out of `cam.r`
  fields (requires the `PT`→`T` conversion — exactly what we want to avoid).
  Decision needed on which the tool supports.
- **O4 — Which `cam.i` fields are honored in a hybrid start?** Verify against
  CESM/ExoCAM source whether `TS1–4`, `TSICE`, `SICTHK`, `SNOWHICE`,
  `ICEFRAC` from `cam.i` are used or overridden by CICE/docn state in an
  active-ice SOM configuration. Determines whether surface-temperature
  acceleration on aquaplanets goes through `cam.i` or is ice/ocean-file-only.
- **O5 — SOM ocean on hybrid.** If hybrid is chosen for atmosphere targets,
  is losing the `somtp` state acceptable (it re-spins from `pop_frc`
  climatology), or must the tool generate a case-specific `pop_frc` /
  fall back to also editing `docn.r`? (Does docn even read `docn.r` on
  hybrid? Verify.)
- **O6 — rpointer & naming discipline.** In-place edits (prototype style,
  with `_backup` copies) keep rpointers valid; writing *new* dated restart
  files would require rewriting `rpointer.*` and case XML. Decide: in-place
  with backups (simple, matches prototype and casemgr's copy-into-rundir
  flow) vs. new-file-set (auditable, no accidental clobber). The
  post-step consistency check (safeguard 6) needs the pre-jump set retained
  either way — rollback depends on it.
- **O7 — Frozen-regime reference set.** Provide a snowball/cold-regime
  restart set before the `aqua_ice` plugin is written (§2, regime note).
- **O8 — cpl.r budget fields.** CESM tracks water/energy budgets in `cpl.r`
  (`budg_*`). A mass-conserving ice jump changes the true inventory the
  budgets describe. Check whether budget diagnostics tolerate the jump or
  need resetting (likely benign — diagnostics only — but verify).

## 6. Recommendation sketch (for discussion, not settled)

Start with the `aqua_ice` target via **in-place `cice.r` (+`docn.r`) edits and
plain continuation**, because (a) it's the prototype's empirically survived
path, (b) no `.i.` files exist for ice anyway, and (c) it defers the
`.i.`-regeneration machinery (O3). Treat **hybrid reboot as the mechanism for
atmosphere-temperature targets** later, after O4/O5 are verified against
source. Either way, the step driver must snapshot the full pre-jump restart
set for the safeguard-6 rollback.
