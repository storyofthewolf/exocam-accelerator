# Jump experiments: what worked, what did not

A running log of every in-model jump and every measurement that changed the method. Add a dated entry for each experiment (newest at the bottom of its section), and update **Lessons so far** when an entry changes what we would do next time. Record negative results as carefully as positive ones. Design detail lives in `docs/ocean-jump.md` and `docs/restart-integration-questions.md`; this file records outcomes.

Entry template: case, restart, what was changed and by how much, the advice it came from, what we expected, what happened (numbers), verdict, lesson.

## Lessons so far

1. **A jump keeps only the part that lies along the slow trajectory.** Sea ice (a slow variable on its own) held almost all of a ×1.5 jump. In the hot runs, the vapor column is the slow reservoir: D4's coupled jump, which filled it at the same time as the ocean, held its whole step (TS flat about +2.7 to +3 K above the pre-jump state through year 80, 19 years). D2's ocean-only probe gave back about 2 K of its first-year +5.3 K within 3 years; at year 80 it reads as an overshoot past D2's equilibrium (see the D2 entry), not a jump that failed to hold, but it is not yet confirmed.
2. **Near-equilibrium parts of the state are overshot by a uniform jump.** The thin substellar ice melted back after ×1.5. Scale only what is still spinning up (the tapered ice jump). In-model (2026-10-07), the tapered jump landed within 1 % of the target thickness in all three conduction-law cases (pt01, pt03 second jumps; pt02 first), and growth afterwards matched the Stefan 1/h rate to within a few percent; the per-cell melt-back has not yet been mapped.
3. **Measure profiles and patterns from time means, never from instantaneous restarts.** Two cam.i snapshots gave the right shape for D5's vertical profile but about twice the amplitude (weather).
4. **The slow vertical response of a hot atmosphere is not its fast response.** In D5 the interannual (fast) response warms the whole column by about 1–2 K per K of surface warming. The decades-long (slow) trend cools the mid-troposphere and warms the upper troposphere about twice as fast as the surface. A jump that imposes the fast pattern will relax back.
5. **Judge convergence from energy_top over several years, never from one year or from energy_bot.** Single-year energy_top in D2 ranges from −18 to +6 W/m²; energy_bot is noisier still.
6. **Advice from a post-jump window is contaminated by the post-jump transient.** The ice runs needed about 20 years to settle; a second jump advised from 18 post-jump years was unreliable.
7. **Interim checks need the run directory.** History files reach the archive only when a job ends, so a mid-run check must build the post-jump years from the run directory and join them to the archived series (only the native columns are valid across the join).
8. **Pick probe cases with low internal variability.** D2 has multi-year TS swings of about ±1 K and single-year energy_top excursions near −18 W/m² (year 68: an eight-month shortwave and longwave event, not a data error). A probe response of a few K cannot be read against that in 10 years; D2 needed 20 years and the user's reading of the Gregory plot (see its entry). Its column cloud water (TGCLDLWP, TGCLDIWP) oscillates with a 5–7 yr period, a candidate source of the swings.
9. **In a hot run, energy_top before a jump includes the vapor column still filling.** D4 had energy_top +3 to +6 W/m² and ran about 6.5 W/m² above energy_bot; after a coupled jump that added 17 % vapor, energy_top fell to about 0 and stayed there for 9 years (TS flat at 373.5 K), and TMQ stopped growing. The advisor's local-curve TS_eq (378 K) overstated the equilibrium; a plain annual N-on-TS line over the last 15–20 pre-jump years (373.1–373.6 ± 0.4 K) got it right.
10. **The coupled cam.r jump is numerically safe in CESM 1.2.1 when TEOUT is recomputed.** D4's first-step energy-fixer correction was −1.4e-4 W/m², with no NaNs, and the model ran on without a test segment. S4 and S5, each with +47 % vapor, also completed their first 10 post-jump years without a crash.
11. **Large-imbalance cases are where jumps pay (user, 2026-10-04).** D4 skipped its slow vapor build-up in one coupled jump. Cases with energy_top above about +5 W/m² (S4, S5) are where the real savings are; an advisor refusal for an unresolved slope is not a reason to wait when the slope is small and the imbalance large (a small λ puts the equilibrium further away, not closer).
12. **An overshoot settles quickly and is a usable result.** S4's coupled +8.7 K overshot its equilibrium, cooled for about 8 years, and has sat at about 337.4 K with energy_top ≈ 0 since (years 81–90 mean +0.5 W/m²): it reached equilibrium about 10 years after the jump, from the warm side. The two-point secant through the pre-jump state and the first post-jump 5-yr mean (334.6 K) undershot that by about 3 K, so it is a bracket, not a prediction.
13. **Near a runaway, the model's numerics become a limit.** S5 warmed about 12 K after its coupled jump with energy_top staying at +8 to +9.5 W/m², then crashed at model month 0082-05 in the dynamical core (MAPZ vertical remap; T max 394 → 1772 K within a few diagnostic prints). A run heading for a much hotter regime may need more conservative dynamics settings to follow it.
14. **Vapor added by a coupled jump scales at about 3.7–4 % per K near 375 K**, consistent with Clausius–Clapeyron (L/(R_v T²) ≈ 3.7 %/K at 375 K) times the profile's mass-weighted gain: D4 +15 % for +3.7 K, D6 +2 % for +0.85 K (gain 0.63).
15. **Some cold cases have no N ≈ 0 equilibrium in this model.** pt04's −2 W/m² is the fusion heat of snowfall accumulating on stagnant thick ice (no seasonal melt, no glacier flow); the climate converged long ago. Check the ice/snow budget before treating a negative N in a partly ice-covered case as spin-up, and keep the snow/ice thickness ratio when jumping ice in such a case (it sits at the flooding threshold). **Revised 2026-10-20:** the TOA imbalance is real and cannot be "corrected" to zero; the planet loses that energy to space for good. Report N_TOA as it is, split by reservoir (atmosphere, slab, ice + snow): pt04 is converged in climate (atmosphere and slab ≈ 0) but not in energy (−2.1 W/m² into the ice reservoir).
16. **An ice jump only removes the part of the deficit that goes through the surface.** Split energy_top into energy_bot (conduction-limited freezing, ∝ 1/hi; a thickness jump removes it) and the gap energy_top − energy_bot (fusion heat of snowfall that never melts; independent of hi). pt01–pt03 (190–203 K, ICEFRAC ≥ 0.94, LHFLX 3–17 W/m²): energy_bot ≈ energy_top, gap 0 to −0.2, and the jumps worked. pt04 (256 K, ICEFRAC 0.68, LHFLX 66): energy_bot −0.2, gap −1.9, and nothing for a jump to remove. The conduction-law offsets a (−0.11, −0.20, −0.61 for pt03, pt01, pt02) rank with evaporation: the same snowfall floor at small amplitude.
17. **Average energy_top over 30+ years in noisy runs.** D1 and D2 (annual σ 3.3 and 4.5 W/m², lag-1 r ≈ 0.3–0.4, a 5–7 yr cloud-water oscillation) show 5- and 10-yr means of +2 to +3 W/m² that are only about 1σ (a 10-yr mean holds about 5 independent years); their 30–40-yr means are +0.4 to +0.7.

## Sea ice (cold aquaplanets)

### 2026-09-28: grp4 pt01 and pt03, cice.r ×1.5 at 0141 — PASS

- **Change:** `vicen` and `eicen` × 1.5 in place (uniform), restart 0141. Advice from the Stefan conduction law `N = a + b/hi`, target removing half the imbalance, capped at 1.5.
- **Expected:** N_after −1.50 (pt01), −1.78 (pt03) W/m².
- **Observed (20 post-jump years):** pt01 landed hi 77.18 m (expected 76.98), N −1.55 against the law's −1.47; pt03 hi 86.54 (86.32), N −1.79 against −1.75. TS within tolerance of the Gregory relation. About 170 years of spin-up skipped (the unjumped runs would reach these thicknesses near year 310).
- **What did not work:** the substellar thin ice (h < 20 m, about 11 % of the area in pt01) was already near local equilibrium; the ×1.5 overshot it, so 5–10 m ice thinned back and 10–20 m ice stalled. N started about 0.2 W/m² more negative than the law and relaxed back over about 20 years.
- **Lesson:** led to the **tapered jump** (2026-09-29): scale each cell by its Stefan growth weight, leave partial-ice cells alone. Retro-tested on the 0141 data (peak 1.5, mean ×1.47, N_after −1.54 against observed −1.55); not yet run in-model.

### 2026-09-29: second-jump advice for pt01/pt03 at 0161 — not jumped

- Advice from the 18 post-jump years fitted `a = +2.70, b = −333` against the pre-jump law's `a = +0.016, b = −117`, which the run still obeys. The post-jump N relaxation contaminates the fit. Not jumped.
- **Open:** wait longer, reuse the pre-jump law, or refuse a positive asymptote.

### 2026-10-04: grp4 all four — pt01/pt03 second jump, pt02 first (tapered, first in-model), pt04 explicit ×1.5 — running

- **State before (10-yr means):** pt01 yr 160 TS 204.5 K, N_TOA −1.5, hi 78.9 m; pt03 yr 160 TS 192.2, −1.7, 88.5 m; pt02 yr 50 TS 201.3, −4.0, 28.7 m (never jumped); pt04 yr 140 TS 254.5, −2.2, 35.1 m (never jumped).
- **pt01/pt03, second jump at 0161, tapered peak ×1.5:** the post-jump-only fit is still contaminated (a = +2.7), so the law was fitted across the first jump on native annual means, years 121–160 (`AdvisorConfig(window_years=40, which="native", detect_jumps=False)`, via the API; no CLI flag yet): pt01 N = −0.212 − 105.6/hi (corr −0.997), pt03 N = −0.194 − 140.8/hi (−0.999); raw factors 2.38 / 2.28, clipped to 1.5. Taper baseline 0151 → 0161 (both after the first jump): pt01 73 % full weight / 12 % unscaled, global hi ×1.464, expected N after −1.11; pt03 76 % / 2 %, ×1.462, N after −1.27. The advice's "years skipped" (~53) is wrong: its h² slope spans the first jump. At the post-jump Stefan rate the step is roughly 250–400 years.
- **pt02, first jump at 0051, tapered peak ×1.5:** `advise` on years 1–50 (int2, 30 yr): N = −0.024 − 115.3/hi, corr −1.000 (the same law as pt01 before its jump); raw 2.01, clipped to 1.5. Taper baseline 0041 → 0051: 70 % full / 15 % unscaled, global hi ×1.468, N after −2.63, skips ~55 yr. pt02 runs about 22 h per 5-yr segment, so this saves about 10 days of compute. **First in-model tapered jump.**
- **pt04, explicit uniform ×1.5 at 0141 (user decision):** the law does not apply. TS (254.5 K) and ICEFRAC (0.684) have been flat for about 100 yr, while hi grows linearly at about 2.1 m/decade and does not slow (corr(N, 1/hi) 0.45). The latent heat of that growth (about 1.4 W/m²) is most of its steady −2.2 W/m². `jump --ice-factor 1.5` (no advice, so `check` cannot score it); hi 35 → 53 m, about 85 yr of buildup at the present rate. Partial-ice cells scaled too.
- Files: advice `$HPC_SCRATCH/grp4_advice/{pt01,pt03}_advice_law121.json`, `pt02_advice.json`, tapers `pt0[123]_taper1004.json` (+ `.taper.nc`); locally `../scratch/grp4/jump1004/`. Backups and jump logs in each `run/exocam_accelerate/`.
- **Run:** in place, RESUBMIT=2 with each case's own STOP_N: pt01 58720425 (yrs 161–190), pt02 58720426 (STOP_N 5: yrs 51–65), pt03 58720427 (161–190), pt04 58720428 (141–170).
- **Watch:** whether the taper avoids the thin-ice melt-back seen after the first ×1.5; whether pt04's growth rate stays unchanged after the step (a thickness-independent sink would show as no change in N).

### 2026-10-07: grp4 after the 2026-10-04 jumps — pt01, pt02, pt03 PASS; pt04 held, imbalance unchanged

- Trends regenerated (`run_trend_batch.sh --cam --cice`) into `$HPC_SCRATCH/grp4_check1007/<pt>`, locally `../scratch/grp4/check1007/` (jump logs alongside); viewer directory `../scratch/grp4/view_1007/`. `check` run locally with each case's `.accel.json`.
- **pt01 (second tapered jump, peak ×1.5, at 0161), years 161–190, 28 settled: PASS.** Landed hi 116.40 m vs 116.28 expected (+0.1 %); N −1.09 vs law −1.10 W/m²; TS 0.32 K below the Gregory relation. hi growth 1.2 → 0.79 m/decade (Stefan 1/h predicts 0.80). At the present rate (d(h²)/dt ≈ 19 m²/yr) 79 → 117 m would take about 400 yr; about 580 yr skipped over both jumps. Years 186–190: TS 203.71 K, N −1.05, hi 118.5 m.
- **pt03 (second tapered jump at 0161), years 161–190, 28 settled: PASS.** Landed hi 130.47 vs 130.32 (+0.1 %); N −1.22 vs law −1.26; TS −0.43 K from the Gregory relation. Growth 1.37 → 0.91 m/decade (Stefan 0.93). About 380 yr skipped. Years 186–190: TS 191.06 K, N −1.20, hi 132.9 m.
- **pt02 (first in-model tapered jump at 0051), years 51–60, 8 settled: PASS.** Landed hi 44.60 vs 44.27 (+0.8 %); N −2.85 vs law −2.55 (diff −0.29, tolerance 0.40). Growth 3.3 → 2.09 m/decade (Stefan 2.1). About 60 yr skipped, as advised. Years 61–65 running.
- **pt04 (explicit uniform ×1.5 at 0141), years 141–170; `check` cannot score it (no law).** Years 111–140 → 143–170: hi growth 2.03 → 1.46 m/decade (Stefan would give 1.29, linear 2.03); TS 254.3 → 255.9 K, a **+1.4 K step at the jump** then flat (unexplained: thicker ice should cool the surface; ICEFRAC 0.684 → 0.679, about 32 % open water); energy_top −2.29 → −2.04, energy_bot −0.38 → −0.13 W/m². Snow depth accelerated: 10.9 m at year 140 growing 0.66 m/decade, to 17.4 m at year 170 (about 1.7 m/decade at first, slowing; snow was not scaled). Hypothesis (untested): pt04's linear ice growth came mostly from snow-ice formation by flooding; thicker ice rides higher, flooding stops, and the snow accumulates instead. The jump held but did not move pt04's steady −2 W/m² sink.
- **Implications:** the conduction-law tapered jump is reliable; the ×1.5 cap limits it (raw factors were 2.0–2.4). pt01's law (asymptote about −0.2 W/m²) puts N = −0.5 at hi ≈ 360 m, thousands of years at present rates, so several more jumps are needed. pt01/pt03 are eligible for a third jump at 0191 (law fitted across the earlier jumps); pt02 follows its pre-jump law, so a second jump could reuse it rather than wait 20 yr for a clean refit; pt04 needs a decision on what its equilibrium is before any further jump. Per-cell maps (thin-ice melt-back under the taper) not yet made.

### 2026-10-08: grp4 pt01 and pt03, third tapered jumps (peak ×1.5) at 0191 — running; pt02 deferred

- **Law fitted across the earlier jumps** (native annual means, years 151–190, `AdvisorConfig(window_years=40, which="native", detect_jumps=False)` via the API; series `$HPC_SCRATCH/grp4_check1007/<pt>/`): pt01 N = −0.204 − 104.1/hi (corr −0.993; before: −0.212 − 105.6/hi); pt03 N = −0.114 − 146.5/hi (corr −0.9986; before: −0.194 − 140.8/hi, so the offset has moved toward zero). Raw factors 2.61 (pt01) and 2.23 (pt03), clipped to 1.5.
- **Taper** (baseline 0181 → 0191, both after the 0161 jump): pt01 72.5 % full weight, 16.2 % partial, 11.3 % unscaled, global hi ×1.471 (vicen ×1.480, eicen ×1.483); pt03 74.8 / 22.9 / 2.3 %, ×1.468 (vicen ×1.479, eicen ×1.482).
- **Change:** pt01 hi 118.7 → ~174.6 m, expected N −1.05 → −0.80 W/m²; pt03 hi 133.1 → ~195.4 m, N −1.21 → −0.87.
- **Years skipped:** the advice says ~70–74, wrong again because its h² slope spans earlier jumps. From d(h²)/dt over years 171–190 (pt01 18.7, pt03 24.3 m²/yr): about **880 yr (pt01) and 840 yr (pt03)**.
- Files `$HPC_SCRATCH/grp4_jump1008/` (`<pt>_advice.json`, `<pt>_taper.json`, `.taper.nc`); pristine backups and jump logs in each `run/exocam_accelerate/`. Pre-flight passed (rpointers at 0191, archived set pristine).
- **Run:** in place, STOP_N=10 nyears, REST_N=10, RESUBMIT=2 (years 191–220): pt01 58780896, pt03 58780897.
- **pt02:** its years 61–65 segment (58758389) was still running; its second jump at 0066 (law across years 26–65, taper 0056 → 0066) is still to do.

### 2026-10-08: pt04 diagnosed — no N ≈ 0 equilibrium in this model: snowfall accumulates on stagnant ice; the ×1.5 jump shifted the climate +1.4 K

> **Corrected 2026-10-20:** the +1.4 K is not from the jump. The unjumped rerun from the same pristine 0141 warmed by the same amount (see the 2026-10-20 pt04 entry). The snowfall budget below stands; the "corrected imbalance ≈ 0" reading does not (lesson 15).

- **Method:** NCO decade means (`ncra`) of cam.h0 and cice.h for years 131–140 (before the jump), 146–155 and 161–170 (`$HPC_SCRATCH/pt04_diag/`, locally `../scratch/grp4/pt04_diag/` with `budget.py`, `budget2.py`). Global means weighted by CICE `tarea` and CAM `gw`; bands by angular distance from the substellar point (lat 2°, lon 160°, the FSDS maximum). CICE namelist: `kdyn = 1` (EVP dynamics on), `kitd = 1`, `tr_pond = .true.`.
- **The TOA deficit is the fusion heat of snow that never melts.** energy_top / energy_bot / gap (W/m²): −2.28 / −0.37 / −1.91 (131–140), −1.98 / −0.10 / −1.88, −2.03 / −0.10 / −1.93. The snowfall onto ice is 216 kg/m²/yr water equivalent (CAM PRECSC+PRECSL and CICE `snow` agree), peaking at 405 mm/yr 70–90° from the substellar point, and its fusion heat is 2.25–2.29 W/m². Snow and ice melt (melts, meltt, meltb, meltl) return only about 0.2 W/m², all at the ice edge (50–90°). So about 2.1 W/m² is sequestered, against the 1.9 W/m² gap: the budget closes within about 10 %. The water is supplied by the open substellar ocean (LHFLX 406 W/m² within 30°), which the slab ocean refills without limit.
- **The ice is effectively stagnant.** Per band, thickness grows at the thermodynamic rate (70–90°: about 0.33 m/yr against `dvidtt` 0.29–0.36), so transport is small; the 60–100 m slabs are too strong to move. (`dvidtd` in the decade means is fill-value garbage; transport was inferred from the thickness change.)
- **Flooding (snow-ice) made most of the ice before the jump:** snoice 0.161 m/yr against congel 0.058 (about 75 %). Snow depth sits at the flooding threshold: hs/hi ≈ 0.31 (0.30–0.31 in every band) against (ρw − ρi)/ρs = (1026 − 917)/330 = 0.33. The ice-only ×1.5 jump dropped the ratio to about 0.20, so flooding stopped: snoice halved to 0.083 (146–155) while snow piled up (hs 10.8 → 15.4 m), then recovered to 0.157 by 161–170 as hs/hi returned to 0.30. The ice-only jump broke the snow/ice ratio; the original `accelerate.py` scaled snow too.
- **The +1.4 K warming is on the night side and comes from the atmosphere,** not the ice edge (ICEFRAC unchanged in every band, so not a partial-ice-cell artifact). ΔTS by band (0–30, 30–50, 50–70, 70–90, 90–120, 120–180°): +0.3, +0.5, +0.6, +1.2, +2.1, +1.8 K; TREFHT the same. Conduction through the thicker ice/snow to the night-side surface fell (fcondtop −0.88 → −0.61 W/m² at 90–120°), which alone would cool it. Dayside cloud liquid thinned slightly (−3 % near the substellar point), absorbed sunlight rose 6 W/m² there (+1.4 W/m² global), OLR +1.1. Why the jump moved pt04 into this state, and whether it is a distinct climate state, is not known; it has persisted 30 years.
- **Conclusion (user agrees, 2026-10-08): a model limitation.** N ≈ 0 needs about 216 kg/m²/yr of snow/ice to melt somewhere; a tidally locked planet has no seasonal melt, and CICE has no glacier-like flow of thick floating ice toward warm water (the sea-glacier mechanism; Goodman & Pierrehumbert 2003). Ice thickness grows about 0.2 m/yr indefinitely. Adding an ice-flow mechanism would be a model change (the user's). pt04's climate had converged before the jump: years 111–140 TS 254.4 K, ICEFRAC 0.684, energy_top −2.3 W/m² with about 2.1 W/m² sequestered as snow/ice (corrected imbalance ≈ 0). If pt04 is continued, roll back to the pristine 0141 first.
- **Method implications (proposals, not built):** for ice cases, report a corrected imbalance energy_top + Lf × (snowfall − melt); refuse a conduction-law ice jump when growth is accumulation-driven (snow-ice dominant, growth independent of h); when an ice jump is made in a snowy case, scale snow with ice (`--snow-factor` equal to the ice factor). S1–S3 (ICEFRAC 0.18–0.66, small negative N, flat TS) probably share this regime; untested.

### 2026-10-08: pt04 re-jump with snow scaled (ice ×1.5 + snow ×1.5) from the pristine 0141 — running (test of the +1.4 K)

- **Question:** was the first jump's +1.4 K step caused by breaking the snow/ice ratio (flooding stopped), or by the added thickness itself? If the former, a jump that keeps hs/hi should not warm; if the latter, it should warm as much or more (snow insulates like ~6–7× its depth of ice).
- **First experiment preserved (user decision: move aside, not overwrite):** years 141–170 history (360 cam.h0, 360 cice.h), restart sets 0151/0161/0171 and the run directory's 0170-12 h0 and 0171 restart set moved into `$HPC_ARCHIVE/exocam_ML_grp4_pt04/jump1_iceonly/{atm/hist,ice/hist,rest,rundir}` (same filesystem, nothing deleted). Its trend series and decade means are in `../scratch/grp4/check1007/pt04/` and `../scratch/grp4/pt04_diag/`.
- **Rollback** (`exocam-accelerate rollback --archive ...`) from the pristine archived 0141 set; the ice-only jump log retired to `.rolledback.json`; rpointers back to 0141.
- **Jump** (`jump --ice-factor 1.5 --snow-factor 1.5`, re-applied from the pristine `.pre-accel.nc`): vicen 171451 → 257176, eicen −5.384e13 → −8.076e13, vsnon 54250 → 81376, esnon −6.574e12 → −9.860e12 (all ×1.500); snow/ice volume ratio kept at 0.316. Partial-ice cells scaled too (as in the first jump), so the comparison differs only in snow.
- **Run:** in place, STOP_N=10 nyears, REST_N=10, RESUBMIT=1 (years 141–160), job 58781369.
- **Compare against the ice-only run** (decade 146–155, and the annual series): global TS (ice-only: +1.42 K over years 131–140), the night-side TS change by band (+2.1 K at 90–120°), snoice (ice-only: halved, 0.161 → 0.083 m/yr), hs/hi (ice-only: 0.31 → ~0.20 → 0.30), energy_top (−2.28 → −1.98), fcondtop and TREFHT by band, and the dayside cloud liquid/absorbed sunlight change.

### 2026-10-16: grp4 checks — pt01/pt03 third jumps PASS (held: close enough); pt04 snow-scaled re-jump warmed +1.4 K too; pt02 second tapered jump at 0066 — running

- Trends regenerated (`run_trend_batch.sh --cam --cice --int1 1 --int2 10`, exocam-trend 12c88ea with per-file caching) into `$HPC_SCRATCH/grp4_check1016/<pt>`, locally `../scratch/grp4/check1016/trend/`; viewer directory `../scratch/grp4/view_1016/`.
- **pt01 (third tapered jump at 0191), years 191–220, 28 settled: PASS.** Landed hi 174.67 m vs 174.59 expected; N −0.76 vs law −0.80 W/m²; TS on the Gregory relation (+0.00 K). Years 216–220 N −0.74.
- **pt03 (third tapered jump at 0191), years 191–220, 28 settled: PASS.** Landed hi 195.47 vs 195.38; N −0.81 vs law −0.87; TS +0.01 K. Years 216–220 N −0.80.
- **User decision:** pt01 and pt03 are close enough to call converged; hold, no fourth jump.
- **pt04 (ice ×1.5 + snow ×1.5 from the pristine 0141), years 141–160:** decade means 131–140 → 146–155: TS 254.57 → 255.96 K (ice-only run: → 255.92), energy_top −2.28 → −1.91 (−1.98), energy_bot −0.37 → −0.01 (−0.10), ICEFRAC 0.68 unchanged, hi 35.1 → 56.0 m, hs 10.8 → 17.2 m. **The warming is the same with the snow/ice ratio kept**, so the +1.4 K is not from stopping flooding; the added thickness itself does it. User: pt04 is problematic — ICEFRAC, qi and hi trends are out of alignment with N and TS. No further action decided.
- **pt02, second jump at 0066, tapered peak ×1.5.** Law fitted across the 0051 jump (native annual means, years 26–65, `AdvisorConfig(window_years=40, which="native", detect_jumps=False)` via the API): N = −0.641 − 99.6/hi (corr −0.995; 30-yr native −0.752 − 95.3/hi, 40-yr int2 −0.630 − 101.5/hi). The offset has moved from the pre-jump law (−0.024 − 115.3/hi), consistent with N sitting ~0.3 W/m² below that law after the first jump. Raw factor 2.88, clipped to 1.5. Taper baseline 0056 → 0066: 69.4 % full weight, 15.7 % partial, 14.9 % unscaled (raw peak 3.19), global hi ×1.474 (vicen ×1.480, eicen ×1.484), hi 47.5 → ~70.0 m; expected N −2.63 → −2.08 W/m². Years skipped: the tool says ~45 (h² slope spans the 0051 jump); from d(h²)/dt over years 56–65 (19.4 m²/yr) about **137 yr**.
- Files `$HPC_SCRATCH/grp4_jump1016/` (`pt02_advice.json`, `pt02_taper.json`, `.taper.nc`); pristine backup and jump log in `run/exocam_accelerate/`. Pre-flight passed.
- **Run:** in place, STOP_N=5 nyears, REST_N=5, RESUBMIT=2 (years 66–80), job 58883658.

### 2026-10-16: pt04 rewound to the pristine 0141 — continuing unjumped for 50 years

> **Corrected 2026-10-20:** the unjumped rerun warmed +1.4 K too, so the jumps did not move pt04 into the warmer state (see the 2026-10-20 pt04 entry).

- **Summary of the two pt04 jumps** (both ×1.5 at 0141, both 20+ years run): ice only (2026-10-04) and ice + snow with hs/hi kept (2026-10-08) gave the same result. Each warmed the planet by +1.4 K, mostly on the night side, made energy_top only about 0.3 W/m² less negative (−2.28 → −1.9 to −2.0) and noisier, and left ICEFRAC at 0.68. The added thickness, not the snow/ice ratio, moves pt04 into a different, warmer state rather than toward its own equilibrium, and N is no closer to zero. pt04's climate had already converged before any jump (years 111–140); its steady −2 W/m² is snowfall/ice sequestered on stagnant ice (lesson 15), which a jump cannot remove.
- **User decision:** the jump added noise and warmed the planet; rewind and continue with no jumps for 50 years and see. No move-aside this time: the jumped run's years 141–160 in the archive (history and restart sets 0151/0161) are overwritten as the rerun reaches them (user: that is fine). Its trend series stay in `../scratch/grp4/check1016/trend/pt04/` and its jump log in `../scratch/grp4/check1016/`. The ice-only run remains in `$HPC_ARCHIVE/exocam_ML_grp4_pt04/jump1_iceonly/`.
- **Rollback** (`rollback --date 0141-01-01-00000 --archive ...`) from the pristine archived 0141 set; rpointers at 0141; the snow-scaled jump log retired to `.accel.rolledback.json`, which replaced the ice-only run's retired log of the same name (a copy is in `../scratch/grp4/check1007/`).
- **Run:** STOP_N=10 nyears, REST_N=10, RESUBMIT=4 (years 141–190), job 58883777.
- **Watch:** whether TS and ICEFRAC stay flat at 254.4 K / 0.684 and energy_top at about −2.3 W/m², with hi still growing about 2 m/decade. Until the rerun has passed year 160, archived years 141–160 mix the jumped and unjumped runs; check file dates before generating trends.

### 2026-10-20: pt02 second jump PASS; third tapered jump (peak ×1.5) at 0081 — running 20 years

- Trends regenerated through year 80 (`trend.py --cam --cice -n 960 --int1 1 --int2 10 --save-data`, in `~/pythonWorkspace/exocam-trend/data/`; locally `../scratch/grp4/check1020/`, viewer directory the same).
- **Second jump (0066), years 66–80, 13 settled: PASS.** Landed hi 70.23 m vs 70.01 expected (+0.3 %); N −1.96 vs law −2.05 W/m²; TS −0.16 K from the Gregory relation; hi +0.14 m/yr. 5-yr means (TS, energy_top, hi): 56–65 199.4 K, −2.74, 46.6 m; 66–70 198.4, −2.06, 70.5; 71–75 198.3, −1.97, 71.2; 76–80 198.2, −1.90, 71.9.
- **Third jump advice** (law across both earlier jumps, native annual means, years 41–80, `AdvisorConfig(window_years=40, which="native", detect_jumps=False)` via the API): N = −0.610 − 99.67/hi (corr −0.994; at 0066 −0.641 − 99.6/hi, so the law holds; 30-yr window −0.396 − 112.3/hi). N now −1.84, raw factor 3.58, clipped to 1.5, uniform N_after −1.53.
- **Taper** (baseline 0071 → 0081): 69.2 % full weight, 15.4 % partial, 15.4 % unscaled (raw peak 4.26), global hi ×1.477 (vicen ×1.483, eicen ×1.486); hi 72.2 → ~107 m, expected N −1.57 W/m², TS ~197.7 K. The tool says ~40 yr skipped (h² slope spans the earlier jumps, so an underestimate).
- Files `$HPC_SCRATCH/grp4_jump1020/` (`pt02_advice.json`, `pt02_taper.json`, `.taper.nc`; local copy of the advice `../scratch/grp4/check1020/pt02_advice_0081.json`); pristine backup and jump log in `run/exocam_accelerate/`. Pre-flight passed; dry run and real jump identical.
- **Run:** in place, STOP_N=5 nyears, **REST_N=1** (was 5: 5-yr segments take 22.3 h of the 24 h limit), RESUBMIT=3 (years 81–100), job 58940519. When done: trends, `check`, and a fourth-jump decision (the law needs hi ~3–4× today's to reach pt01/pt03's −0.8 W/m²).

### 2026-10-20: pt04 — the +1.4 K is the model's own; closed as converged in climate; energy budget and snow enthalpy checked

- **Unjumped rerun from the pristine 0141, years 141–190 (job 58883777 and resubmits), trends `../scratch/grp4/check1020/`:**

  | years | TS (K) | energy_top (W/m²) | ICEFRAC | hi (m) |
  |---|---|---|---|---|
  | 131–140 | 254.57 | −2.28 | 0.684 | 35.1 |
  | 141–150 | **256.03** | −2.02 | 0.680 | 37.1 |
  | 161–170 | 256.02 | −2.15 | 0.678 | 41.0 |
  | 181–190 | 255.89 | −2.03 | 0.679 | 44.8 |

  hi continues smoothly from 35 m, so the archive holds only the unjumped run. TS rose by the same +1.4 K as both jumped runs (255.9–256.0 K), ramping over months 141-07 to 142-02 (no step at the restart), then flat for 50 years. **The jumps did not cause the warming; they had essentially no effect on TS.** The original run stopped at 140, so there was never an unjumped continuation to compare against until now. Namelist and setup sections of atm, cpl, ice and ocn logs from the last original segment (131–140, 2026-09-27) and the rerun (141–150) differ only in restart names, dates and step counts. Not checked: whether the executable changed between July and September.
- **User decision: pt04 closed.** Converged in climate at about 256.0 K (the post-141 state, held 50 yr, longer than the 30 yr of the 254.6 K state). The two quasi-stable states about 1.4 K apart are not explained; if pt04's equilibrium is published, compare decade-mean maps of TS, cloud, ICEFRAC and hs for 131–140 vs 141–150 first.
- **Why pt04 has no N ≈ 0 state and pt01–pt03 do** (last-decade means, W/m² unless stated):

  | | pt03 | pt01 | pt02 | pt04 |
  |---|---|---|---|---|
  | TS (K) | 190 | 203 | 198 | 256 |
  | ICEFRAC | 1.00 | 0.985 | 0.945 | 0.68 |
  | LHFLX | 3.4 | 9.9 | 17 | 66 |
  | TMQ (kg/m²) | 0.08 | 0.23 | 0.49 | 14.3 |
  | hs/hi | 0.007 | 0.009 | 0.016 | 0.31 |
  | energy_top | −0.82 | −0.75 | −1.94 | −2.03 |
  | energy_bot | −0.84 | −0.63 | −1.73 | −0.16 |
  | gap | +0.02 | −0.12 | −0.21 | −1.87 |

  See lesson 16. The gap is the fusion heat released in the atmosphere when vapor freezes to snow and radiated to space; the ocean gave the water up as vapor and gets it back as ice, so the ice + snow system carries that debt until the snow melts, which on pt04 it never does.
- **Energy diagnostic (user, 2026-10-20):** no "corrected imbalance" holds. N_TOA = dH_atm/dt + dH_slab/dt + dH_ice+snow/dt is the planet's real energy loss. Report N_TOA with that split; climate convergence = TS flat and the atmosphere and slab terms ≈ 0; energy equilibrium = N_TOA ≈ 0. pt04 meets the first, not the second. Proposed tool change (not built): report an ice-reservoir uptake term beside N_TOA, never subtracted from it.
- **Snow-enthalpy check (one-time):**
  - Code (pt04 SourceMods, read only): CICE `ice_therm_vertical.F90` adds new snow with `qsnew = −rhos·Lfresh` (lines ~3919–3929), and the column energy check (~4465) includes `−fsnow·Lfresh` and stops the model on failure. docn SOM (`docn_comp_mod.F90` 646–656) charges `−(snow + ioff)·latice` for snow on open water. Slab hblt 50 m everywhere, qdp 0.
  - Global closure, years 141–150 vs 181–190 (`ncra` decade means of cice.h in `$HPC_SCRATCH/pt04_enthalpy/`, locally `../scratch/grp4/pt04_enthalpy/` with the docn restarts 0141/0191 and the SOM forcing file): grid-mean hi 37.12 → 44.81 m, hs 11.40 → 13.94 m, so the fusion debt is −Lf·(ρi dhi + ρs dhs)/dt = **−2.09 W/m²**, plus an estimated −0.15 to −0.25 for the new ice and snow being colder than 0 °C (Tsfc about −24 °C global mean); slab +0.03 (mean somtp 277.87 → 278.10 K); atmosphere ≈ 0 (TMQ 14.39 → 14.33). Total −2.2 to −2.3 against N_TOA −2.12 (years 141–190): closes within 5–10 %, inside the uncertainty of the sensible estimate. A missing snow enthalpy would show as a ~2 W/m² mismatch. CICE's `qi`/`qs` history fields (units "J") do not behave as column enthalpy (`qs` falls while hs grows 22 %) and were not used.
- **Model improvements discussed (out of scope for this tool; ExoCAM CICE SourceMods, the user's call):**
  1. Cap ice + snow thickness and calve the excess, with its enthalpy, to open-ocean cells near the ice edge, where the slab must melt it (as CLM caps snow on land ice and sends frozen runoff to the ocean). Closes energy and water exactly; the cap is arbitrary, the melt location is not neutral and may move the ice edge.
  2. Sea-glacier flow: gravity spreading of thick floating ice (Weertman, velocity ~h³) toward thinner ice, melting at the edge (Goodman & Pierrehumbert 2003; Pollard & Kasting 2005). Physical steady thickness and melt location; research-grade (top thickness category, enthalpy/snow conservation, substellar geometry).
  3. Ocean heat transport under the ice (a prescribed day-to-night q-flux in the SOM; input data, no code change) melts from below; it changes the climate on its own and does not guarantee closure.
  Also: track the water inventory if night-side water trapping becomes a question. In this model it cannot happen: the slab is an unlimited water source, and CICE ice thickness is unrelated to the 50 m slab depth.

## atlasfu aquaplanets (D-series, 4 bar, 340–377 K, ice-free; S-series, 265–372 K, S1–S3 partly ice-covered)

### 2026-10-01: D2, ocean-only probe, docn.r somtp +8 K at 0061 — PASS at year 80: overshot to the cooling side, continuing

- **Change:** somtp +8.000 K uniform in place (351.885 → 359.885 K), probe mode, probe-dt 4, assumed ocean heat fraction 0.5, so a TS step of +4 K. Advice `$HPC_SCRATCH/atlasfu_d2_probe/D2_probe4_hr2.json`. Continuation job 58680001, years 61–70; archived through restart 0071.
- **Final check (2026-10-02, years 61–70, 8 settled):** verdict WAIT.
  - Landed: first-year TS +5.33 K net of drift against +4.00 expected, so the implied ocean heat fraction was 0.67, not 0.5.
  - TS: 357.87 (61), 356.78, 354.92, 353.79, 354.01, 354.24, 355.42, 355.61, 355.38, 353.80 K (70). Pre-jump drift +0.163 ± 0.054 K/yr; post-jump (years 63–70) +0.07 ± 0.22.
  - Against the pre-jump trend extrapolated (years 41–60 fit), D2 averages +0.8 K ahead over years 64–70 (about 5 years of spin-up), but year 70 is 0.5 K *behind*. With multi-year swings of about ±1 K, the net gain is between zero and about 1 K: not distinguishable from zero with confidence.
  - energy_top +3.4 → −2.2 ± 3.3 W/m², λ = 2.4 ± 1.4 W/m²/K: restoring, no sign of runaway, not resolved. Year 68 alone has energy_top −18.1 W/m²; without it the post-jump mean is +0.05.
  - Year 68 is a real event, not a data error: energy_top −5 to −30 W/m² in every month, first from higher OLR (FLNT +15 to +31 W/m² in Mar–Jul), then from less absorbed sunlight (FSNT −16 to −44 W/m² in Aug–Nov), with TS about +1 K throughout. TS then fell 1.6 K in year 70. D2 had a similar 1.5 K dip in year 58 before the jump.
  - The atmosphere followed the ocean: in year 61 the whole column warmed (+15 K at 144 hPa, +3 K at 660 hPa, +10 K at 2850 hPa; 0.5–2.7 K per K of TS), and over years 64–70 it stays 3–5 K warmer than years 56–60 at every level, with water vapor up 10–95 %. D2 (TS ≈ 355 K) is in the uniform-warming regime that D5 left near TS ≈ 372 K.
- **Interpretation (not proven):** the jump is restored by the fast feedback (λ ≈ 2.4 W/m²/K, relaxation within about 3 years), while the run's slow drift, with N ≈ +3.5 W/m² sustained for decades and a flat Gregory relation, continues at its own rate. The fast λ restoring a perturbation is not the slow effective λ along the trajectory. A state jump only helps if it moves whatever sets the slow drift, which we have not yet identified for D2.
- **Lesson (year 70; revised at year 80, below):** the ocean-only probe was safe (no crash, no runaway) but bought little. D2 is a poor probe target because of its internal variability.
- **Decision (2026-10-02, user):** extend D2 by 10 years in place (job 58691769, years 71–80) to decide whether ocean-only somtp jumps are worth doing at all. Rollback not indicated.
- **Year-80 data (2026-10-03):** trends regenerated through year 80 with TMQ, TGCLDLWP, TGCLDIWP (`$HPC_SCRATCH/atlasfu_d2_probe/trend80w`, locally `../scratch/atlasfu/d2check80w`). Global means from `runmgr check --energy`: years 76–80 TS 353.7 K, energy_top −2.8 W/m²; years 73–80 TS 353.8 K, energy_top −0.9 W/m². exocam-trend 10-yr means at year 80: energy_top −0.40, energy_bot +1.04 W/m².
- **Check at year 80 (2026-10-03, years 61–80, 18 settled):** `check` verdict PASS (probe branch: equilibrium bracketed). TS 352.31 → 354.08 K; energy_top +3.38 → −1.27 ± 2.91 W/m², λ 2.6 ± 1.6 W/m²/K. TS drift +0.163 ± 0.057 K/yr before, −0.114 ± 0.088 after (+2.7σ change); relaxation fit τ 11 ± 4 yr, TS_eq 352.8 ± 0.7 K. TMQ jumped with the probe (≈1300 → 1450–1500 kg/m²) and has not relaxed. TGCLDLWP oscillates 55–155 g/m² with a 5–7 yr period (TGCLDIWP in phase), present before and after the jump; a candidate source of D2's multi-year TS and energy_top swings (untested).
- **Reading (user, 2026-10-03):** the probe jumped D2 to the cooling side of the Gregory relation, i.e. past equilibrium, and the run is cooling into it; post-jump years in the Gregory plot and every time series agree. The agent had read the same numbers as "little gained"; the user's reading is that the probe bracketed the equilibrium margin. Confirm or refute with several more decades: refuted if the long-term subtle warming trend resumes.
- **Viewer (2026-10-03):** the Gregory plot now colours post-jump years on their own red scale (years after the jump), unconnected, with a post-jump int2 line. Every post-jump year sits on the cooling side (energy_top < 0 at TS above the pre-jump state), moving back toward the pre-jump Gregory line; TS, energy_top and energy_bot time series agree.
- **Revised lesson (year 80):** an ocean-only probe large enough to cross the equilibrium is informative: the run's response from the far side (cooling at −0.11 K/yr toward TS_eq ≈ 352.8 K) brackets the equilibrium from above, which the pre-jump drift (+0.16 K/yr, energy_top +3.4) alone did not show. Whether it saved spin-up time depends on that equilibrium being right; the agent's "falls back" reading at the time missed this.
- **Decision (2026-10-03, user):** continue D2 in place, STOP_N=10, RESUBMIT=4 (years 81–130), job 58708736. Refuted if the slow warming trend resumes.

### 2026-10-02: D5 vertical warming profile, measured before the coupled jump — no jump yet

- **Measurement:** per-level global means of T on model levels (exocam-trend `--profile T,Q`), annual means, years 1–60. Gain = per-level warming per K of surface warming.
- **First attempt (cam.i 0051 → 0061):** mid-troposphere −6 K/K, model top +20 K/K. Rejected as weather-contaminated; replaced by time means.
- **Time means (trend over 10 and 20 yr, two annual means over 10 and 20 yr):** all agree. Upper troposphere (100–400 hPa) +1.3 to +2.3 K/K; zero crossing near 600 hPa; mid-troposphere (800–3500 hPa) cooling, strongest at about 2000 hPa (−2.6 ± 0.5 K/K over years 41–60); lowest two levels +0.5 to +0.9.
- **Robustness (2026-10-02):**
  - Per-level trends over years 41–60 have t-statistics of 4–8 at both ends of the profile, with AR(1)-inflated errors.
  - Not a coordinate artifact: model levels drift to higher pressure as vapor adds mass (up to +4 hPa/yr near the surface), but the gain at fixed pressure (log-p interpolation) is the same within 0.2 K/K.
  - The pattern is a **regime change**: in 10-yr windows ending at years 20, 30 and 40 the column warmed nearly uniformly (gain about 0.5–1.5 everywhere). The mid-troposphere cooling appears from about year 40 (TS ≈ 372 K) and persists through year 60. In absolute terms, between years 31–35 and 56–60 the 2300 hPa level cooled 6.9 K while TS rose 5.8 K and 100 hPa warmed 11.9 K.
  - Fast and slow responses differ: regressing detrended interannual anomalies on TS anomalies gives +1 to +2 K/K at every level; the slow trend is the cooling.
  - Q rises at every level (+5 to +30 % per decade).
  - Physical reading (hypothesis, not tested): increasing upper-level water-vapor absorption of stellar radiation warms the upper troposphere and starves the mid-troposphere of solar heating.
- **Consequence:** the current `atm-profile` cuts the jump off where warming first turns negative (about 4000 hPa), leaving a mass-weighted gain of 0.06. The coupled jump as built is effectively a small ocean-only jump. Following the measured profile would require cooling part of the column, which the code does not do. **Open decision.**
- Outputs: `$HPC_SCRATCH/atlasfu_d5_jump/prof60/`, locally `../scratch/atlasfu/d5view/prof60/` and `trend60p/`.
- Untested: whether D4 or D6 show the same regime change.
- **Decision (2026-10-02, user):** D5 not jumped; continued unjumped for 5 × 10 yr (job 58691775, RESUBMIT=4) alongside D1 and D3 (58691773, 58691774) to spin up naturally. D4 untouched; D6 continues. A coupled ocean+atmosphere jump remains to be tried on a case still in the uniform-warming regime.

### 2026-10-03: D4, first in-model coupled jump (docn.r somtp + cam.r T/q) at 0061 — preliminary: holds, near balance at ~373.5 K (PASS at year 70, see 2026-10-04)

- **Advice:** `advise-ocean` local curve on energy_top, last 40 yr: TS 370.78 K, energy_top +4.05 ± 0.49 W/m², α_diff 0.63 ± 0.20 W/m²/K (resolved), TS_eq ≈ 378.3 K, ocean heat fraction 0.29 (ocean-only advice: somtp +13.1 K, clipped to +10). `atm-profile --method trend --window 20` (years 41–60): gain ≈ 1 K/K at the surface, max 1.78, mass-weighted 0.96, se ≤ 0.10 below 200 hPa; trend-10 and annual-10 agree (0.88, 0.73). D4 is still in the uniform-warming regime (no mid-troposphere cooling as in D5). Files `$HPC_SCRATCH/atlasfu_jump1003/D4_{ocean,coupled}.json`; constants from the archived `atm.log.260930-161123`.
- **Change:** TS +3.75 K (370.78 → 374.54); somtp +3.753 K uniform (370.60 → 374.35); cam.r T +3.69 K at the lowest level (max +6.67), q at fixed RH: +616 kg/m² vapor (+15 %), PS +60.3 hPa, +72.6 W yr/m²; TEOUT recomputed (constants reproduce it to 1.2e-4). Skips ~19 yr. Expected energy_top after ≈ +2.4 W/m².
- **Run (user decision 2026-10-03):** the 1-month test (job 58709230) was cancelled before it started; D4 runs 10 yr straight from the jumped restart instead (job 58709261, STOP_N=10 nyears, RESUBMIT=0), with REST_OPTION nmonths/1 left on (monthly restart sets, ~97 MB each). "If it crashes it crashes." Check atm.log (energy fixer first step, NaNs, PS) afterwards or on a crash; restore REST_OPTION nyears/10 for later segments.
- **Model health (job 58709261, 8.4 h in, model 0065-09):** first-step energy-fixer correction −1.4e-4 W/m² (later steps 1e-3–1e-2), no NaNs in atm.log, global-mean PS 4474.6 hPa at the first step, 4480.5 hPa at 0065-08. About 0.56 model yr per wall-clock hour.
- **Interim result (run directory, years 61–65; area-weighted annual global means from cam.h0):**

  | year | TS (K) | energy_top (W/m²) | TMQ (kg/m²) | PS (hPa) |
  |---|---|---|---|---|
  | 56 | 370.33 | +5.90 | 4067 | 4398.6 |
  | 57 | 370.25 | +4.51 | 4122 | 4403.9 |
  | 58 | 370.45 | +3.85 | 4134 | 4405.1 |
  | 59 | 370.62 | +4.01 | 4161 | 4407.7 |
  | 60 | 370.78 | +2.90 | 4183 | 4409.9 |
  | 61 | 373.81 | −1.26 | 4887 | 4478.9 |
  | 62 | 373.47 | +0.13 | 4894 | 4479.5 |
  | 63 | 373.52 | −0.25 | 4894 | 4479.5 |
  | 64 | 373.59 | −0.12 | 4888 | 4479.0 |
  | 65 (8 mo) | 373.41 | −0.27 | 4902 | 4480.4 |

  - Landed: TS +3.0 K in year 61 against +3.75 planned, then flat near 373.5 K for 4.7 years (pre-jump warming about +0.1–0.2 K/yr).
  - energy_top fell from +3–6 W/m² to about 0 (−0.3 to +0.1 after year 61), against an expected +2.4. TMQ rose 17 % (dry run: 15 %) and has not grown since; before the jump it grew about 25 kg/m² per year.
- **Interpretation (preliminary, not proven):** D4 may be at or near equilibrium at about 373.5 K, below the advisor's TS_eq ≈ 378 K. Before the jump energy_top ran about 6.5 W/m² above energy_bot (heat going into the atmosphere, mostly the vapor column); with the column filled, that gap and energy_top both closed. If it holds, the coupled jump skipped D4's slow vapor build-up, more than the ~19 years the advice predicted. Five years is short against the ±1 K multi-year swings seen in D2; the rest of the segment and `check` (energy_bot, implied heat ratio) will tell.
- **Next:** when job 58709261 ends, regenerate D4 trends through year 70 and run `check`; restore REST_OPTION nyears/10 (and realign restarts to January) before further segments; purge or thin the ~120 monthly restart sets.
- Same pre-analysis (2026-10-03), not jumped: S5 (363.5 K, +8.5 W/m², α_diff 0.59 ± 0.33 unresolved; coupled probe +5.7 K = +49 % vapor; 10-yr profile window cut by a noisy mid-troposphere dip, 20-yr gives 1.61) and S4 (331.6 K, +6.2 W/m², α_diff 0.20 ± 0.31, last-10-yr λ −0.43 ± 0.39 in the 310–330 K cloud-transition regime; coupled probe +8.7 K = +47 % vapor).

### 2026-10-04: D4 at year 70 — PASS, at equilibrium; D2 years 81–90

- **D4 (coupled jump at 0061), job 58709261 completed years 61–70.** `runmgr check --energy`: years 66–70 TS 373.5 K, energy_top +0.5 W/m²; years 62–70 TS 373.5 K, energy_top +0.3 W/m². exocam-trend (`~/pythonWorkspace/exocam-trend/data/exocam_atlasfu_D4_0001-01-0070-12_*`): settled years 63–70 energy_top +0.11 W/m², TS drift −0.003 K/yr. **Reading (user): a picture of success — a jump to N = 0 that holds.** The jump skipped D4's remaining spin-up.
- **Check:** the first check of the year-70 data said FAIL (N 2.03 W/m² below the pre-jump Gregory line, tolerance 0.75): it scored the run against a line that overestimated the equilibrium. Changed (5ea283b): energy_top within noise of zero and TS flat over the settled years is PASS ("reached equilibrium"), whatever the pre-jump line said. D4 now reads PASS.
- **OPEN (user, 2026-10-04): revisit the verdict criteria as a whole** — PASS / FAIL / WAIT for checks and jump / no jump / probe for advice. The present gates (local-curve α_diff resolution, line correlation ≥ 0.9, linear-vs-saturating agreement, the Gregory-line tolerance in `check`) disagree with the user's reading of the Gregory plots and int1 series on D4, S4 and S5. `--override-gate` does not bypass the line-quality filters.
- **Monthly restarts:** no action needed. After a successful run CESM keeps the last restart set and cleans up the rest (user, 2026-10-04); the earlier note to purge them was wrong. Restore REST_OPTION nyears/10 before D4's next segment.
- **D2 (ocean probe at 0061), years 81–90 (job 58708736; years 91–100 running as 58716887):** years 81–90 TS 353.4 K, energy_top +0.7 W/m²; years 86–90 TS 354.0 K, energy_top +1.0 W/m². Above the year-80 relaxation TS_eq (352.8 ± 0.7 K); not yet distinguishable from the 5–7 yr cloud oscillation. Trends in `~/pythonWorkspace/exocam-trend/data/exocam_atlasfu_D2_0001-01-0090-12_*`.
- Viewer directory with D2, D4, S4, S5: `../scratch/atlasfu/view_1004/`.


### 2026-10-04: S4 and S5, coupled jumps (docn.r somtp + cam.r T/q) at 0061 — running

- **Decision (user, 2026-10-04):** jump both. Cases with energy_top above +5 W/m² are where jumps save the most, and the user reads linear Gregory lines on both. The advisor's full step refused both (local-curve α_diff unresolved: S4 0.20 ± 0.31, S5 0.59 ± 0.33), and its line fit refused S4 (correlation and linear-vs-saturating gates), so both jumps are probe-sized: TS step = recent trend × years, coupled (ocean heat fraction 1).
- **Basis (annual means, years 56–60):** S4 TS 330.4 K, energy_top +6.2 W/m², TS trend +0.53 K/yr; annual N-on-TS line over 20–30 yr: λ 0.20–0.24 W/m²/K, TS_eq 354–359 K (± 8–10). S5 TS 362.7 K, energy_top +8.5 W/m², trend +0.45 K/yr; λ 0.23–0.33, TS_eq 388–400 K (± 13 or more). Neither jump overshoots unless λ is above about 0.7 (S4) or 1.7 (S5).
- **S4:** `S4_coupled.json` (probe 15 yr × 0.58 K/yr = +8.70 K; profile `--method trend --window 20`, years 41–60, mass-weighted gain 0.84). somtp 328.63 → 337.33 K; cam.r T +8.91 K at the lowest level (max +8.99); vapor +198 kg/m² (+47 %); PS +19.4 hPa; +33.7 W yr/m². TEOUT reproduced to 7.2e-4.
- **S5:** `S5_coupled55.json` (`advise-ocean --probe --probe-dt 5.5`, then `atm-profile --method trend --window 20`, years 41–60, mass-weighted gain 1.61). somtp 363.10 → 368.60 K; cam.r T +5.64 K at the lowest level (max +11.57); vapor +1341 kg/m² (+47 %); PS +131.4 hPa; +160.1 W yr/m². TEOUT reproduced to 1.4e-4. Skips ~14 yr by the trend.
- Advice and profiles in `$HPC_SCRATCH/atlasfu_jump1003/`; constants from the archived `atm.log` copies there (`atm.log.S4`, `atm.log.S5`); pristine backups and jump logs in each `run/exocam_accelerate/`.
- **Run:** in place, plain continuation, STOP_N=10 nyears, REST_OPTION nyears/10, RESUBMIT=4 (years 61–110): S4 job 58718608, S5 job 58718609.
- **Expect (by analogy with D4):** energy_top falls well below its pre-jump value as the filled vapor column stops absorbing heat. Check first: atm.log first-step energy-fixer correction, NaNs, PS.

### 2026-10-04: S1, S2, S3 at year 60 — near equilibrium, not jumped; continued

- **Trajectory (annual means, years 51–60):** S1 TS 265.2 K, energy_top −1.2 W/m², ICEFRAC 0.66, hi 12.5 m; S2 274.5 K, −0.7, 0.57, 8.4 m; S3 285.0 K, −0.3, 0.18, 3.3 m. TS has been flat for 20 years (trends −0.08, −0.04, +0.06 K/decade, none significant); ICEFRAC steady (S2 −0.003/decade).
- **Ice thickness is still growing** steadily: S1 +1.1 m/decade, slowing as the Stefan law expects; S2 +0.74 m/decade, not slowing; S3 +0.25. The latent heat this releases (ρ_i L dh/dt × ICEFRAC ≈ 0.7, 0.4, 0.05 W/m²) matches energy_bot and most of the small negative energy_top, so the run is not still cooling. Thickness is the one slow variable left; it barely affects TS or N.
- **Decision (user):** no jump. Continued in place, STOP_N=10 nyears, RESUBMIT=3 (years 61–100): S1 58720301, S2 58720302, S3 58720303. Trends in `~/pythonWorkspace/exocam-trend/data/exocam_atlasfu_S[123]_0001-01-0060-12_{cam,cice}.txt`. D4 relaunched by the user (job 58718868).

### 2026-10-05: status of the D- and S-series — S4 overshot, S5 shows no restoring response, D4 energy_top drifting up

All 11 atlasfu cases running (jobs: D1 58722355, D2 58725083, D3 58722362, D4 58726640, D5 58722505, D6 58722340, S1–S3 58720301–3, S4 58726514, S5 58726619). Numbers are area-weighted global means from `runmgr check --energy` over the archived output (the running segment is not included); "derived" values are the 10-yr mean minus the last-5-yr mean, not measured directly.

- **S4 (coupled +8.70 K at 0061; pre-jump years 56–60 330.4 K, energy_top +6.2 W/m²), years 61–70 (job 58718608): overshot.**

  | window | TS (K) | energy_top (W/m²) |
  |---|---|---|
  | years 61–65 (derived) | ~340.5 | ~−0.8 |
  | years 66–70 | 338.3 | −5.4 |
  | years 68–70 | 337.7 | −3.2 |
  | year 70 | 337.6 | −4.6 |
  | years 61–70 | 339.4 | −3.1 |

  energy_top went strongly negative and TS is falling about 0.5 K/yr, so the equilibrium lies below about 338 K, far below the pre-jump line's TS_eq of 354–359 K. The jump log said it would not overshoot unless λ exceeded about 0.7 W/m²/K; the apparent slope from pre-jump to years 66–70 (ΔN ≈ −11.6 W/m² over ΔTS ≈ +7.9 K) is about 1.5. As in D4 (lesson 9), part of the pre-jump +6.2 W/m² was probably the vapor column filling, so the pre-jump N-on-TS line was not λ. S4 started in the 310–330 K cloud-transition regime, where λ was negative over the last pre-jump decade, so the slope may be state-dependent too (untested). TS is still about 7 K above the pre-jump state, approaching its equilibrium from the warm side.
- **S5 (coupled +5.50 K at 0061; pre-jump 362.7 K, +8.5 W/m²), years 61–70 (job 58718609): no restoring response.**

  | window | TS (K) | energy_top (W/m²) |
  |---|---|---|
  | years 61–65 (derived) | ~370.3 | ~+8.4 |
  | years 66–70 | 371.7 | +9.0 |
  | years 69–70 | 371.9 | +9.3 |
  | year 70 | 372.0 | +9.5 |
  | years 61–70 | 371.0 | +8.7 |

  TS rose about 9 K in total (jump plus continued warming at about 0.3–0.4 K/yr) while energy_top stayed at +8.4 to +9.5 W/m² and is creeping up. Unlike D4 and S4, energy_top did not fall when the vapor column was filled. Under the runaway rule this is a warning ("runaway greenhouse suspected"), verdict WAIT; rollback is the user's call. D4 and D5 are balanced at 374–377 K, so S5's +9 W/m² at 372 K differs from the D-series (S5's configuration not compared here).
- **D4 (coupled jump at 0061; PASS at year 70), years 71–80 (job 58718868):** years 71–75 (derived) ~373.6 K, ~+0.2 W/m²; years 76–80 373.8 K, +1.4; years 78–80 373.9 K, +1.8; year 80 373.9 K, +2.0. TS is flat but energy_top has drifted up by about 1–2 W/m² over the last 3–5 years. This could be the 5–7 yr cloud oscillation seen in D2; watch the next segment (years 81–90, job 58726640).
- **Unjumped, 5-yr means (TS, energy_top):** D1 years 86–90 338.8 K, −0.8 W/m²; D2 (ocean probe at 0061) years 96–100 353.5 K, +0.8; D3 years 86–90 356.9 K, −0.1; D5 years 86–90 376.7 K, +0.7 (converged on its own after the 2026-10-02 no-jump decision); D6 years 51–55 374.3 K, **+2.6** (the only D case still clearly out of balance; behind the others); S1–S3 years 56–60 unchanged (265.2 / 274.5 / 285.1 K; −1.1 / −0.6 / −0.1), first segment of years 61–100 still running.
- **Annual series (exocam-trend, native columns; `$HPC_SCRATCH/atlasfu_check1005/{S4,S5,D4}`, locally `../scratch/atlasfu/view_1005/`):**
  - S4 TS (K) / energy_top (W/m²), years 60–70: 331.60/+7.44, **340.55/+3.46**, 341.03/+0.37, 340.65/−5.35, 340.04/−2.94, 339.79/−0.85, 339.94/−8.76, 338.72/−8.84, 337.81/−4.21, 337.61/−1.14, 337.56/−4.79. Years 61–65 mean 340.41 K, −1.06; 66–70 338.33 K, −5.55. TS peaked in year 62 and its fall is slowing (−1.2, −0.9, −0.2, −0.0 K/yr over years 67–70) while energy_top stays negative. Two-point secant between years 56–60 (330.44 K, +6.23) and 66–70 puts N = 0 near 334.6 K (slope 1.49 W/m²/K).
  - S5 years 60–70: 363.50/+9.61, **369.80/+7.84**, 370.28/+7.20, 370.39/+8.26, 370.67/+9.54, 370.90/+8.20, 371.36/+8.23, 371.43/+8.87, 371.68/+8.67, 371.85/+8.85, 372.04/+9.28. energy_top dipped only about 1 W/m² in year 61 (D4: about 4) and is back to the pre-jump level. TS warming slowed from about 0.45 K/yr (years 51–60) to about 0.25 K/yr (years 61–70) at the same imbalance, consistent with a larger column heat capacity after the vapor increase; not yet a runaway signature in the rate, but no restoring response in N.
  - D4 years 71–80: energy_top +0.23, −0.37, +0.03, −0.73, +0.69, +0.02, **+1.38, +1.64, +1.20, +1.87**; TS 373.45–373.60 through year 76, then 373.73–374.00. Years 71–75 mean 373.53 K, −0.03. Since year 77, TS and energy_top have risen together, which is not a restoring response; possibly the cloud oscillation seen in D2 (TGCLDLWP not yet checked for D4).
- **Next:** look at these in the viewer (Gregory plots), then the user's decision on S5 (continue or roll back to 0061) and whether S4 is left to settle.

### 2026-10-07: S5 crashed, S4 settled, D4 converged; D1/D3 stopped, D5 continued, D6 coupled jump

Means from `runmgr check --energy` (archived output) and, for D1/D3/D5/D6, exocam-trend series through the last archived year (`$HPC_SCRATCH/atlasfu_trend1007/`, locally `../scratch/atlasfu/trend1007/`, with `summarize.py`; viewer `../scratch/atlasfu/view_1007/`). Uncertainties are AR(1)-inflated.

- **S5 (coupled +5.5 K at 0061): FAILED at model month 0082-05** (2026-10-06 17:42, segment years 81–90 of job 58718609's resubmits). `MAPZ_MODULE` abort (`mapz error - k0found ... -1`: vertical remap failed, layer interfaces crossed); global T max had been 394–408 K for the whole segment, then 698, 666, 1772 K over the last three diagnostic prints; U max 228 m/s. Years 71–80: TS 373.8 K, energy_top +9.4 W/m²; years 76–80 374.5 K, +9.5. About 12 K warmer than before the jump, with no restoring response in energy_top over 20 years. The crash came about 21 years after the jump, in a run still warming; the hot regime is the more likely cause than the jump edit itself, but the two cannot be separated. Pristine restart for any continuation: archived 0081-01-01.
  - **Reading (user, 2026-10-07):** D4 and S5 are both successful jumps. S5's jump raised the temperature without reducing the imbalance; the signs point to S5 being in a runaway, or at least transitioning towards a much hotter regime. Continuation (with more conservative dynamics) or not: open.
- **S4 (coupled +8.7 K at 0061): settled at about 337.4 K.** Years 81–90 TS 337.4 K, energy_top +0.5 W/m²; years 86–90 337.4 K, −0.9; year 90 337.5 K, −0.3. TS flat since about year 69. Running years 91–100 (58764896).
- **D4 (coupled jump at 0061): converged.** Years 91–100 373.6 K, +0.2; 96–100 373.5 K, −0.1. The upward drift of years 77–80 was transient. Next segment queued (58766834).
- **D2:** years 111–120 354.1 K, +0.7; 116–120 354.6 K, −3.0 (111–115 about +4.4): large swings continue, 10-yr mean near balance. Running years 121–130, the last resubmit (58764590).
- **S1–S3:** years 71–80 265.2 / 274.5 / 285.2 K; −0.9 / −0.6 / −0.2 W/m²; running years 81–90.
- **D3: converged; stopped (COMPLETE at year 110).** TS 356.8–357.1 K for 50 yr; 30-yr trend +0.05 ± 0.02 K/decade; energy_top 30-yr mean +0.11 ± 0.12; TMQ flat. Equilibrium 356.9 K.
- **D1: converged; stopped (COMPLETE at year 110).** TS 338.4 ± 0.4 K for 50 yr (30-yr trend −0.18 ± 0.22 K/decade); energy_top 30-yr mean +0.35 ± 0.58. Gregory lines over 15–40 yr: TS_eq 338.3–338.6 K. Years 106–110 +3.1 W/m² is one swing of a large internal oscillation (5-yr means −2.7 to +3.1; cloud liquid water roughly anti-correlated), not a trend. Equilibrium about 338.4 K; quote its balance as a 30–50-yr mean.
- **D5: slow drift; continued unjumped 2 × 10 yr (job 58775076, years 111–130, REST_OPTION nyears/10).** TS +0.2 K/decade over 30 yr; energy_top 20-yr mean +1.0 ± 0.2; TMQ +50 kg/m²/decade. Budget closes: vapor storage (L·dTMQ/dt ≈ 0.4 W/m²) + energy_bot (≈ 0.6) ≈ energy_top. Gregory 15–30 yr: TS_eq 377.4–378.0 K (weak, r −0.1 to −0.4). `advise-ocean` refuses (no restoring local curve). Not jumped: the remaining 0.5–1 K is small, and D5's mid-troposphere-cooling profile (2026-10-02) cannot be followed by the coupled jump.
- **D6: coupled full-step jump at 0065-01-01, in place; continuing (job 58775077, years 65–94, STOP_N 10, RESUBMIT 2; REST_N was 1, now 10; STOP_N was 9, now 10).**
  - Before (years 60–64): TS 374.79 K, energy_top +1.42 W/m², energy_bot +0.64, TMQ 4011 kg/m² growing 95 kg/m²/decade; energy_top decaying smoothly (+23 at years 15–19, +2.6 at 45–49).
  - Advice (`advise-ocean --n-fraction 1`, local curve on energy_top, 40 yr): α_diff 1.47 ± 0.13 W/m²/K, TS at N = 0 375.74 K, ocean heat fraction 0.34, relaxation time ~17 yr. Plain annual Gregory lines agree: 30 yr λ 1.48 (r −0.93), 40 yr 1.56 (r −0.99), TS_eq 375.6–376.1 K for 15–40-yr windows. Full step chosen (user's call that half steps are too conservative); the default half step was +0.42 K.
  - Profile (`atm-profile --method trend --window 20`, years 45–64; the 10-yr window is refused as noise): gain 0.97 K/K at the lowest level, 1.3–1.6 above 500 hPa, a dip to 0.3–0.5 ± 0.3–0.5 at 1300–3500 hPa (possibly the onset of D5's regime, still positive); mass-weighted 0.63.
  - Change: somtp 374.818 → 375.668 K (+0.850 uniform); cam.r T +0.82 K at the lowest level (max +1.24); vapor +97 kg/m² (+2 %); PS +9.5 hPa; +11.0 W yr/m². TEOUT reproduced to 1.1e-4 (constants from the archived `atm.log.261004-184048`). Dry run and real jump identical. Files `$HPC_SCRATCH/atlasfu_jump1007/` (`D6_ocean.json`, `D6_coupled.json`, `.atmprofile.nc`, `atm.log.D6`); backups and jump logs in D6's `run/exocam_accelerate/`.
  - Expect energy_top ≈ 0 at TS ≈ 375.7 K (or lower, as in D4, if filling the column lowers N further). Skips roughly 35–50 model years. When the segment ends: atm.log first-step energy fixer, NaNs, PS; then `check`.


### 2026-10-15: atlasfu status — ten of eleven cases converged or closed; D6 restarted after a wall-clock kill

10-yr means from `runmgr check --prefix exocam_a --energy -n 10` (archived output). D1 and D2 also from exocam-trend series (`$HPC_SCRATCH/atlasfu_trend1015/D2`, D2 without per-level profiles; D1 reuses `atlasfu_trend1007/D1`; locally `../scratch/atlasfu/view_1015/` with `summarize.py`). Uncertainties are AR(1)-inflated.

| case | last 10 yr | TS (K) | energy_top (W/m²) | verdict |
|---|---|---|---|---|
| D1 | 101–110 | 338.0 | +1.9 | converged (stopped at 110) |
| D2 | 121–130 | 354.1 | +2.4 | converged within its noise; stopped (COMPLETE at 130) |
| D3 | 101–110 | 356.9 | +0.2 | converged (stopped at 110) |
| D4 | 111–120 | 373.5 | +0.2 | converged; coupled jump at 0061 held about 60 yr |
| D5 | 121–130 | 377.3 | +0.4 | converged without a jump |
| D6 | 55–64 | 374.7 | +1.8 | WALLCLOCK after the 0065 jump; restarted by the user |
| S1 | 91–100 | 265.1 | −0.9 | converged |
| S2 | 91–100 | 274.4 | −0.6 | converged |
| S3 | 91–100 | 285.2 | −0.2 | converged |
| S4 | 101–110 | 337.0 | −0.4 | converged; coupled jump at 0061 overshot, then settled |
| S5 | 71–80 | 373.8 | +9.4 | suspected runaway (user); crashed at 0082-05, not continued |

- **D1:** TS 5-yr means 337.8–338.8 K for 50 yr. 30-yr trend −0.18 ± 0.22 K/decade, energy_top 30-yr mean +0.35 ± 0.58. Gregory lines over 15–40 yr give TS_eq 338.3–338.6 K (r ≈ −0.5). The 10-yr +1.75 rides on one swing (years 96–100 −2.7, 106–110 +3.1). Equilibrium about 338.4 ± 0.4 K.
- **D2:** TS 5-yr means 353.2–354.6 K since year 86; energy_top 5-yr means −3.1 to +4.2. Mean energy_top +2.24 ± 0.40 over 10 yr, +1.40 ± 1.05 over 20, +0.71 ± 0.90 over 30. dTS/dt over 30 yr +0.33 ± 0.41 K/decade; TMQ +31 ± 29 kg/m²/decade, a slow vapor storage like D5's. Gregory lines over 15–40 yr give TS_eq 354.1–354.6 K (weak lines, r ≈ −0.5, but the windows agree). Equilibrium about 354.4 ± 0.5 K. Not jumped: the remaining ≤ 0.5 K is below D2's ±1 K interannual swings and could not be verified. Quote its balance as the 30-yr mean. Narrowing that to about ±0.5 W/m² would take about 30 more years.
- **S1, S2:** persistent small negative energy_top with TS unchanged since year 60. Probably the same snowfall sequestration on stagnant ice as pt04 (lesson 15); not checked.
- **S5:** the user reads it as a runaway (or a transition to a much hotter regime). It is recorded as such and not continued. A continuation from 0081-01-01 with more conservative dynamics would be a separate experiment, worthwhile only if the runaway trajectory is itself a result.
- **D6:** the first post-jump segment (job 58775077, STOP_N 10, REST_N 10) ended WALLCLOCK on 2026-10-10 with nothing archived after year 64. D6 had previously run 9-yr segments. The cause was not diagnosed; the likely reading is that 10 yr did not fit the wall-clock limit, and with REST_N = STOP_N no restart had been written. The user restarted it. Next: atm.log first-step energy fixer, NaNs, PS, then `check` against TS ≈ 375.7 K.
- **Operational lesson:** keep REST_N well below STOP_N (5 or 1) on segments near the wall-clock limit, so a timeout keeps its progress.
- **Tooling note:** exocam-trend with `--profile T,Q` over 130 yr exceeded the 10-min login-node watchdog (killed at 575 s); without profiles it took 3.5 min. Long profile runs belong in a batch job.

### 2026-10-20: atlasfu — D1/D2 swings are variability; D6 jump held, 20 more years running; S5 continued and still running away

10-yr and 5-yr means from `runmgr check --prefix exocam_atlasfu --energy -n 10` / `-n 5` (archived output).

| case | window | TS (K) | energy_top 10-yr / 5-yr (W/m²) | status |
|---|---|---|---|---|
| D1 | 101–110 | 338.0 | +1.9 / +3.2 | converged (30-yr mean +0.35) |
| D2 | 121–130 | 354.1 | +2.4 / +2.2 | converged (30-yr mean +0.71) |
| D3 | 101–110 | 356.9 | +0.2 / +0.0 | converged |
| D4 | 111–120 | 373.5 | +0.2 / +0.1 | converged |
| D5 | 121–130 | 377.3 | +0.4 / +0.5 | converged |
| D6 | 78–87 | 375.4 | +0.6 / +0.9 | jump held; running years 88–107 |
| S1–S3 | 91–100 | 265.1 / 274.4 / 285.2 | −0.9 / −0.6 / −0.2 | converged |
| S4 | 101–110 | 337.0 | −0.4 / −0.3 | converged |
| S5 | 86–95 | 377.5 | +10.4 / +10.0 | continued from 0081; still warming |

- **D1 and D2 are converged; their +2 to +3 W/m² 5- and 10-yr means are variability (user agrees after the viewer, `../scratch/atlasfu/view_1015/`).** Annual energy_top over the last 40 yr: σ 3.3 (D1) and 4.5 (D2) W/m², lag-1 r 0.37 and 0.31, so a 10-yr mean has about 5 independent years and a 1σ of about ±1.5 (D1) and ±2 (D2). Means over the last 5/10/20/30/40 yr: D1 +3.1/+1.75/+0.4/+0.35/+0.4; D2 +2.1/+2.2/+1.4/+0.7/+0.6. 5-yr blocks swing both ways (D1 −2.7 to +3.1; D2 −3.1 to +4.2) with σ about 20 % above white noise. Annual energy_top is anti-correlated with TGCLDLWP (r −0.52 D1, −0.48 D2): every strongly negative block is a high-cloud-water block (D1 96–100: −2.7, 103 g/m²; D2 106–110: −2.4, 105; 116–120: −3.1, 104). TS 5-yr means stay flat (D1 337.8–338.8, D2 353.2–354.6 K). D2 still edges up about +0.3 K/decade (20-yr mean +1.4 above the 30-yr +0.7); a residual drift of a few tenths of a K cannot be excluded but is inside its ±1 K swings. Lesson 17.
- **D6 (coupled +0.85 K jump at 0065): held, near balance at 375.4 K.** Trends through year 87 (`trend.py --cam -n 1044`, locally `../scratch/atlasfu/view_1020/` with both jump logs). 5-yr means (TS, energy_top): 60–64 374.79, +1.42 (before); 65–69 375.54, +0.45; 70–74 375.51, −0.01; 75–79 375.44, +0.02; 80–84 375.32, +0.18; 83–87 375.38, +0.70. Years 66–87: 375.44 K, +0.26 ± 0.10 W/m² (annual σ 0.5, lag-1 r −0.05; far quieter than D1/D2); TS −0.10 K/decade, TMQ −16 kg/m²/decade, a small settling below the jump point. No overshoot: energy_top went from +1.4 to ~0 at once and stayed within −0.6 to +1.5 in every year; years 85 and 87 (+1.0, +1.5) look like D4's transient 77–80 swing. Target was N = 0 near 375.7 K; skipped roughly 35–50 model years. Not scored with `check`. Two failed jobs since the restart (58877159 after 54 s on 10-15, 58893839 after 4 h on 10-17) were not diagnosed; later 5-yr segments with REST_N 1 ran cleanly (12.4 h each).
  - **User decision:** run 20 more years to let it settle. `runmgr continue`, STOP_N=5 nyears, REST_N=1, RESUBMIT=3 (years 88–107), job 58940170.
- **S5: continued by the user from 0081-01-01** (dynamics settings unchanged from D6: nsplit 32, nspltrac 4, nspltvrm 4, div24del2flag 42, ATM_NCPL 48; the earlier values were not recorded). Running past the 0082-05 crash point without trouble: job 58920002 at 0098-09 on 2026-10-19, one resubmit left. TS / energy_top: years 56–60 362.7 K / +8.5 (before the jump); 76–80 374.5 / +9.5; 86–95 377.5 / +10.4; 91–95 378.1 / +10.0. TS is still rising about 0.25 K/yr and energy_top has not fallen at all; at 378 K S5 is +10 W/m² out of balance where D5 balances at 377.3 K. The runaway signature holds.
- **Open:** D6 after year 107 (trends, `check`, close); S5, how far to follow it.

### 2026-10-22: S5 turns over at ~380 K — cloud brightening against falling OLR; 10 natural years before any second jump

- **S5 finished years 96–105** (job 58920002, COMPLETE 2026-10-21, restart 0106-01-01). `runmgr check --energy`: 101–105 380.2 K / +7.0 W/m²; 96–105 379.6 / +8.2. Only one jump ever (coupled at 0061, 2026-10-04; `run/exocam_accelerate/` holds just the 0061 cam.r/docn.r logs and backups); the 0081 restart continued the jumped trajectory unedited.
- **Annual series** (`trend.py --cam -n 1260 --save-data`, locally `../scratch/atlasfu/view_1022/` with the 0061 jump logs). TS / energy_top: 86 376.36 / 11.37; 90 377.50 / 11.32; 95 378.49 / 10.81; 97 378.92 / 9.92; 98 379.17 / 6.96; 100 379.61 / 9.65; 102 380.07 / 6.40; 103 380.10 / 6.17; 104 380.32 / 6.61; 105 380.56 / 5.10. Years 102–105 are the four lowest since the jump. TS still +0.25 K/yr.
- **The turnover is a cloud-shortwave feedback racing a falling OLR, not D1/D2-style variability.** 5-yr means, FSNT / FLNT / TGCLDLWP: 61–65 323.8 / 315.6 / 153; 81–85 321.2 / 309.9 / 174; 91–95 316.6 / 306.8 / 214; 101–105 307.1 / 300.3 / 262. OLR keeps falling as TS rises (~3 W/m² per 5 yr; the runaway signature); absorbed SW now falls faster (−5.2 per 5 yr) as liquid cloud water rises steadily (+70 % since year 85; in D1/D2 it fluctuated about a flat mean). The local N-vs-TS slope went from ~0 (to year 97) to ~2–3 W/m²/K; a line through the last years gives TS_eq ~384–386 K, poorly constrained because the curve is still bending.
- **Year-61 jump in hindsight:** it advanced TS ~5.5 K (~15–20 yr of warming at the pre-jump rate); energy_top could not fall because there was no restoring feedback at 363–378 K.
- **Decision (user):** wait for natural years before a second jump — a jump transient would mask whether the cloud brightening keeps outpacing the OLR decline, and the slope is changing too fast to size a jump. `runmgr continue --set STOP_N=5 --set RESUBMIT=1` from 0106 (REST_N unchanged), job 58967787 → years 106–115. Then: regenerate trends; if N keeps falling, a lower-bound coupled jump of ~1.5–2 K (cannot overshoot unless λ > ~3); if N returns to ~+9, S5 is back to the runaway reading.
