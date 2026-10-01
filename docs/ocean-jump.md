# Hot, ice-free runs: the som_ocean (docn.r somtp) jump

Session 2026-09-30. Target regime: slab-ocean aquaplanets with thick, hot
atmospheres (atlasfu D1–D6: 4 bar CO2, TS 340–375 K, no sea ice) that take many
decades to converge. Offline drivers (outside the repo):
`../scratch/atlasfu/ocean_hindcast.py`, `../scratch/atlasfu/ocean_endtoend.py`.
Trend data: exocam-trend `--cam`, 50 yr (D5 40, D6 10).

## What docn.r holds

From ExoCAM `docn_comp_mod.F90` (SOM mode):

- one variable, `somtp(gsize)`, the prognostic mixed-layer temperature, **in K**;
- `gsize = ni*nj` of the docn domain, in MCT global-index order (lon fastest),
  so `somtp.reshape(nj, ni)` is the lat-lon map. Verified on D4: that reshape
  is smooth in both directions (the transposed one is not), and its area mean
  (369.2 K) and 10-yr rate (0.315 K/yr) match the trend TS (369.1 K, 0.33 K/yr);
- on restart the first coupling step copies `somtp` straight into `So_t`, so an
  edited somtp is the surface temperature the atmosphere sees from step one;
- the docn domain file is named in the run directory's `docn_ocn_in`
  (`domainfile = ...`); atlasfu uses `domain.ocn.4x5_100120_aquaplanet.nc`
  (46×72, all ocean).

`exocam-accelerate somtp-map` writes `somtp(lat, lon)` as an ordinary netCDF.

## Finding 1: the equilibrium is energy_bot = 0; the TOA gap is atmospheric storage

The TOA imbalance (`energy_top`, what `runmgr check --energy` reports as Etop)
exceeds the surface imbalance (`energy_bot`, net flux into the slab). In the
atlasfu runs the gap is ~1–2 W/m² in D1–D3 and 6–9 W/m² in D4/D5 over their
last 15 years.

*Corrected 2026-09-30 (same day).* I first read the gap as an energy leak,
because `energy_top − energy_bot` did not regress cleanly on dTS/dt. That
regression cannot separate storage whose capacity grows with T from a leak
that grows with T. An **equilibrated 3-bar CO2 run** (`scratch/CO2_2949120ppm.dat`,
59 yr, 303 → 375 K) settles it:

- the gap is 30–38 W/m² while TS rises fast and closes to 0.36 W/m² once TS is
  flat (years 35–59): it is heat the atmosphere stores while it warms, not a
  leak;
- the atmosphere's heat capacity (gap / dTS/dt, or the slope of the cumulative
  `energy_top − energy_bot` against TS) grows from ~4 at 345 K to ~16–21
  W·yr/m²/K at 369–374 K. A Clausius–Clapeyron estimate of the vapor column's
  latent heat at 375 K gives ~15, plus ~1 for 3 bar of dry CO2: it is the
  water-vapor reservoir;
- using that C_atm(T), the storage it implies explains the atlasfu gaps to
  ±0.7 W/m² (D1–D3, late D5); D4 has 2–4 W/m² more, plausibly its larger 4-bar
  atmosphere.

Consequences:

- `energy_bot` is still the right coordinate: the slab's equilibrium is where
  the net surface flux vanishes, and the storage term C_atm(T)·dTS/dt bends a
  TOA Gregory line. In the 3-bar hindcast the `energy_bot` line puts the
  equilibrium within +0.5 K of the true 375.1 K from year 16 (+1–2 K at years
  12–14); the `energy_top` line is refused until year 28. On atlasfu the TOA
  line put it 1.5–7 K too hot.
- **An ocean-only jump lands at C_ocean / (C_ocean + C_atm) of its size.** The
  advisor now measures `heat_ratio = C_total / C_ocean` from the cumulative
  heat uptake (∑energy_top, ∑energy_bot against TS) and scales the somtp
  increment by it, clipped to [1, 4]. Measured: ~1.4 (D1, 338 K), ~1.5 (D2),
  ~2.2 (D3, 356 K), ~3.4 (D4, D5 and the 3-bar run at 369–375 K). At 370 K a
  somtp increment of 10 K buys ~3 K of TS; the extra heat goes to the vapor
  column within weeks, a transient SST overshoot of the jump size.
- The mixed layer shows up as C_ocean ≈ 6.4–7.4 W·yr/m²/K (≈ 50 m, `hblt`).

`runmgr check --energy`'s Etop therefore lags convergence (it includes the
storage term) but is not biased at equilibrium.

## Finding 2: the TS(N) relation steepens in the hottest runs

D4/D5 warm steadily (0.33 and 0.55 K/yr over the last 15 years) while
`energy_bot` falls only slowly. (The equilibrated 3-bar run shows the same
steepening near its end — λ_surf falls from 1.5 to 0.8 W/m²/K between 360
and 375 K — and still converges smoothly: steepening is not a runaway.) The local slope dTS/dN_bot steepens from ~−1.3
(20–40 yr windows) to −2.7…−4.9 K per W/m² (last 10–20 yr): sensitivity rising
with temperature. A long window's line then passes *below* the current state,
and naive advice for D4 is a −0.4 K (cooling) jump in a run that is still
warming. The advisor therefore requires the current state (mean of the last
5 native years) to lie on the fitted line within max(0.5, 2σ/√5) W/m²; it
tries windows longest first and refuses when none qualifies. D4 and D5 are
refused on this ground as of year 50 / 40 — correctly, by the settled design
("refuse near a feedback threshold"). Between-period estimates put D4's
equilibrium near ~376 K, well above the long-window line (368 K).

Where the line is accepted but curves, jumps err on the conservative side: at
arrival at the recommended TS, `energy_bot` is still +2 W/m² above the
predicted N_after in D4/D5 (the run warms further than the line says).

## Finding 3: D1–D3 are within noise of equilibrium

At year 50 the recommended jumps are +0.21 K (D1), −0.40 K (D2), −0.18 K (D3):
all below the interannual TS scatter, which the advisor now uses as the floor
for a jump worth making. D2 swings ±2 K from year to year. D6 (10 yr) is too
short to fit.

## The jump

- Coordinate `energy_bot`; linear Gregory line `TS = c0 + c1·N` over the
  longest trusted window (the saturating form alongside; disagreement = kink,
  refused); `c1 < 0` required.
- Target: remove `--n-fraction` (0.5) of the current imbalance;
  `dTS = c1·(N_target − N_now)`; somtp increment `= dTS × heat_ratio`,
  clipped at `--max-dt` (10 K; hard bound 25 K per cell).
- The somtp increment is `dTS × heat_ratio`, the ratio measured over the fit
  window (Finding 1; `--heat-ratio` overrides; the window is reported as
  `heat.window`). The reported "now" imbalance (`N_now`) is the mean of the
  last `recent_years` (5) native years with its standard error (`N_now_se`);
  the last single year is kept as `N_last_year` (on D2 it read -0.08 against a
  noisy 5-yr mean of -0.75 ± 4). `check` reports the ratio the
  first post-jump year implies, which tests the measurement.
- Uniform by default. `pattern` scales each cell by its warming rate between
  two archived docn.r (3×3 box-smoothed, weights clipped to [0, 3], area mean
  1). On D4 (0041→0051) the weights span 0.54–1.56, std 0.14: the warming is
  close to uniform.
- Refusals: ice anywhere in the latest year (early icy years are just skipped),
  no stable relation, current state off the line, jump smaller than the
  interannual scatter.
- `years_skipped = τ·ln(N_now/N_after)` with τ = −c1·C_eff (one box).

Mechanics as for the ice jump: in-place edit of the run-directory copy of
docn.r (`rpointer.ocn` names it), pristine backup and JSON log in
`run/exocam_accelerate/`, the same pre-flight, `check`, `rollback` and
`restore`. Near-freezing cells (possible ice) and masked cells are never
changed.

## Probe mode (decided with the user, 2026-09-30)

The user disagreed with refusing D4/D5: flattening of the imbalance does not
imply a runaway, and a warm perturbation maps the trajectory far more cheaply
than running hundreds more years. The data agree that the recent relation is
simply unconstrained: over the last 10 yr D4's `energy_bot` changes by
−0.1 ± 0.9 W/m² while TS rises 3 K (λ = 0.03–0.10 ± 0.5), D5's λ is 0.2–1.0.

`advise-ocean --probe` therefore offers a jump that is *not* sized to an
equilibrium: TS is stepped ahead by the run's own recent trend ×
`--probe-years` (15) — a TS step; somtp moves by the measured heat ratio times
that, the ratio measured over the longest of 40/30/20/10 yr with usable data
(not the 15-yr trend window: on D2 the ratio reads 1.02 over 15 yr and 2.0
over 40 yr, matching the TOA-surface gap / dTS/dt) — i.e. forward Euler behind the settled time-domain
trustworthiness gate and the hard clip (or `--probe-dt` explicitly, which
skips the gate — D2's ±2 K swings fail it). It reports the smallest λ the
check will resolve after 5 settled years and the λ above which the probe
overshoots (`N_now/ΔT`). With the heat ratio (~3.4) the 10 K clip binds at year 50: D4
+10 K somtp for +2.9 K of TS, D5 +10 K for +3.0 K; a +3 K probe on D2 only resolves λ ≳ 3
(σ of energy_bot ≈ 7 W/m²).

The probe check (`check.check_ocean_probe`) compares the settled post-probe
years with the 5 pre-probe years — the lever arm the natural run lacks:
λ = −ΔN/ΔTS with its standard error; FAIL if the imbalance *grew* in the
direction of the probe (no restoring feedback: the runaway signature — roll
back) or the probe never landed; PASS once λ is significant, with
TS_eq = TS_post + N_post/λ and the side (still short of equilibrium, or
overshot = bracketed); WAIT while the response is inside the noise (PASS with
an upper bound on |λ| after 12 settled years). The viewer's probe toggle draws
the outcome fan: where the post-probe state lands for each λ.

### TS-trajectory readout of the probe check (2026-10-01)

On a noisy run energy_bot is a poor gauge: D2 (353 K) has an annual σ of
7.6 W/m², so after 5 settled years the lever-arm λ resolves only
λ ≳ 2.3 W/m²/K while D2's λ is probably < 1. Annual TS scatters by only
~1 K, and its drift (0.16 ± 0.04 K/yr before the probe) is measurable on both
sides of the probe. `check_ocean_probe` therefore also reads the probe from
TS (`check.ts_relaxation_readout`), as *additional* metrics (`ts_*`); the
energy_bot readout is unchanged.

One-box relaxation, `dTS/dt = (TS_eq − TS)/τ`: two (rate, temperature)
points on one line,

    r_pre  at T_pre   (OLS over the last 20 pre-probe years; at least 8)
    r_post at T_post  (OLS over the settled post-probe years, last 15 at most)
    τ = (T_post − T_pre)/(r_pre − r_post),   TS_eq = T_pre + r_pre τ,
    λ = C/τ  (C = advice `C_eff_W_yr_m2_K`, else `heat.C_ocean`; the
              ProbeCheckConfig.heat_capacity overrides)

Method choice: a straight-line fit over a window short against τ gives the
derivative and the temperature at the window's mid-time (bias ~(L/τ)²/24,
1.5 % for L = 10 yr, τ = 25 yr), so the post-probe transient is handled
without assuming its shape, and (T, r) pairs have Gaussian errors that
propagate by the delta method. A joint exponential fit (TS_eq, τ, landed step)
to a nearly straight 20-yr pre window is badly conditioned and its errors are
not Gaussian; it uses the same information (the change in rate across a known
change in temperature).

Errors: slope se = σ_eff/√Sxx, mean se = σ_eff/√n, with σ the pooled residual
scatter of the two fits (post-window curvature only inflates it) and
σ_eff = σ √((1+ρ)/(1−ρ)), ρ the lag-1 autocorrelation of the pre-probe
residuals clipped to [0, 0.5] (a 20-yr ρ is itself noisy; the cap stops one
lucky run inflating the error without bound). τ and TS_eq errors are first
order (delta method) and are reported only once the rate change is resolved
(≥ 2σ); until then only a lower bound on τ (and the matching upper bound on
λ). The uncertainty of C is not propagated. The pre window starts no earlier
than the advice's `since_year`.

Readout states and the verdict rules (conservative; FAIL rules unchanged
plus one):

- FAIL: probe did not land; non-finite data; energy_bot rose in the probe's
  direction by ≥ 2σ (runaway); **new:** the TS drift *rose* after the probe by
  ≥ 3σ (r_post > r_pre: TS accelerating away — runaway signature). The TS
  rule uses 3σ because a false alarm costs a rollback.
- PASS: at least one readout is significant — energy_bot (≥ 2σ), or TS
  (rate change ≥ 2σ, which yields τ, TS_eq ± se, λ and the side; or TS
  falling back against the probe, r_post < 0 at ≥ 2σ: equilibrium below the
  probe level, bracketed; r_post still ≥ 0 with a resolved rate change:
  equilibrium above) — and no FAIL rule fires and the two do not
  contradict.
- Contradiction: both readable, putting the equilibrium on different sides of
  the probe level (energy "hotter still" vs TS "bracketed", or the reverse):
  WAIT, with the reason printed. Keep running.
- WAIT: fewer than 5 settled years, neither significant, or a contradiction.
  After 12 settled years with neither readable: PASS with bounds
  (|λ| below what energy resolves, τ above the TS lower bound).

Resolution (synthetic D2-like probe: +4 K from 352.7 K, TS_eq 358 K, τ 25 yr,
C = 7, annual σ_TS 1.0–1.2 K, σ_N 7.6 W/m²; 100 noise realisations, readable
fraction by settled years): τ = 25 yr: 16 % at 5, 42 % at 8, 60 % at 10, ~80–95 %
at 12–15; τ = 15 yr: 30 % at 5, 80 % at 8, 92 % at 10; τ = 40 yr: 11 % at 5,
31 % at 10, ~60 % at 15 (rest report a τ lower bound). AR(1) noise with
ρ = 0.3 slows this by 1–3 years. The energy_bot readout alone does not
resolve these cases at all inside 12 years. The pulls of τ and TS_eq are
centred (|median z| < 0.6) in `tests/test_probe_ts.py`.

## Post-jump check (`check`, som_ocean logs)

- landed: first post-jump annual TS moved by the expected
  `somtp_dT/heat_ratio`, net of the pre-jump drift and of the one-box
  relaxation during that year; FAIL if < 30 % and the shortfall exceeds
  2σ of interannual TS;
- on-line: settled-year `energy_bot` against the Gregory line at the observed
  TS, relative to the pre-jump offset; tolerance max(0.75 W/m², 2σ/√n);
- non-finite post-jump data is a FAIL, never a PASS.

## The coupled jump: ocean + atmosphere (decided 2026-09-30)

A somtp-only jump buys little at 365–375 K — the atmosphere drains 2/3–3/4
of it into the vapor column — and deeper, hotter atmospheres will be worse.
The user decided: jump the atmosphere as well, through **cam.r in place +
continuation**, with a **measured vertical profile** and **q at fixed RH**
(`atmos.py`, `restart.apply_atmos_jump`, `cli atm-profile`).

What cam.r needs (CAM FV / ExoCAM cesm1.2.1, verified on atlasfu D4 0051):

- `PT` is T_v/pkz; with the run's own constants (read from `atm.log`: CPDAIR,
  RAIR, ZVIR, SURFACE GRAVITY, PTOP) T is recovered exactly. D4 is
  3.8 bar N₂ + 0.2 bar CO₂ (cp 1024.26, R 288.70, κ 0.2819, zvir 0.5986),
  with 3855 kg/m² of precipitable water;
- `Q` and `DELP` are moist. Adding vapor adds mass (CAM's `dme_adjust` does
  the same): each layer keeps its dry mass, `DELP' = DELP(1−q)/(1−q')`, and
  `PS` rises. Condensate (`CLDLIQ`, `CLDICE`, `LCWAT`) keeps its mass. On a
  restart the dry-mass fixer only prints (`dryairm`, `nlres`), so nothing
  undoes this;
- **`TEOUT` must be recomputed.** The global energy fixer
  (`check_energy_gmean`) applies (energy now − TEOUT) as a uniform heating on
  the first step; with a stale TEOUT it would remove the whole jump at once.
  CAM's formula (dry static energy cp·T + g·z + Φs, kinetic, (L_v+L_i)·q,
  L_i·liq) reproduces D4's stored TEOUT to 0.005 % rms; the jump refuses when
  the constants reproduce it worse than 0.1 % (wrong constants for the case);
- `TCWAT`, `QCWAT` (previous-step T, q of the stratiform scheme) and
  `T_TTEND` shift with T and q.

The profile (`atm-profile`): horizontal-mean T per level from two archived
cam.i (the run writes one with every restart set), per K of area-mean somtp
warming between them, 3-point smoothed, clipped to [0, 4], applied from the
surface up to where the measured warming first turns negative, tapered over
two levels above. D4 (0041→0051, ΔTS 3.16 K): 1.0 K/K at the lowest level,
up to 2.8 K/K near 50 hPa, stratosphere cooling (−4.7 K/K at the top, not
applied). Measured q rises 3.8 %/K near the surface — the fixed-RH
Clausius–Clapeyron rate at 372 K.

The coupled advice keeps the ocean advice's TS target but drops the heat
ratio: somtp moves by ΔTS itself and the atmosphere by gain × ΔTS. Consistency
check: a 1 K coupled jump on D4 puts 22–29 W·yr/m²/K into the atmosphere
(mostly vapor: 203–245 kg/m² per K); the trend data imply C_atm ≈ 18 — the
same order; the first real jump's `check` (implied heat ratio ≈ 1 if the
atmosphere is right) will tell which is closer.

Scale of a jump: the default D4 probe (+4.9 K surface) would add 1328 kg/m²
of vapor (+34 %), raise PS by 130 hPa (+3 %) and the upper troposphere by
14 K. The jump refuses a vapor increase above 50 %; a first real coupled
jump should be smaller (e.g. `advise-ocean --probe-dt 1.5`).

Rehearsed end to end on copies of D4's real restarts in a mock run directory
(`advise-ocean --probe` → `atm-profile` → `jump` → `rollback`): both files
written and verified, rollback restores both bit-identically and retires both
logs. A jump that edits two files is one jump: `check` scores the docn.r log,
`rollback` restores and retires the whole group.

## Not done / open

- No real ocean jump has been made yet; the first one is the Tier-1/2 test.
  D4/D5 are the probe candidates (the Gregory advisor refuses them).
- H2-dominated envelopes: the plugin assumes CAM's constants describe the
  gas (read from atm.log) and liquid-water saturation above 273 K; check
  both before using it on a new atmosphere.
- No coupled jump has been run in the model yet. First real test: a small
  coupled probe on a clone (Tier-1 null jump at zero size first), watching
  atm.log for the energy fixer's first-step heating (`nstep, te` lines: the
  heating should stay at its usual size) and the first months of TS/PS.
