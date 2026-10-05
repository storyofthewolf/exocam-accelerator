# Jump experiments: what worked, what did not

A running log of every in-model jump and every measurement that changed the method. Add a dated entry for each experiment (newest at the bottom of its section), and update **Lessons so far** when an entry changes what we would do next time. Record negative results as carefully as positive ones. Design detail lives in `docs/ocean-jump.md` and `docs/restart-integration-questions.md`; this file records outcomes.

Entry template: case, restart, what was changed and by how much, the advice it came from, what we expected, what happened (numbers), verdict, lesson.

## Lessons so far

1. **A jump keeps only the part that lies along the slow trajectory.** Sea ice (a slow variable on its own) held almost all of a ×1.5 jump. In the hot runs, the vapor column is the slow reservoir: D4's coupled jump, which filled it at the same time as the ocean, held its whole step (TS flat about +2.7 to +3 K above the pre-jump state through year 80, 19 years). D2's ocean-only probe gave back about 2 K of its first-year +5.3 K within 3 years; at year 80 it reads as an overshoot past D2's equilibrium (see the D2 entry), not a jump that failed to hold, but it is not yet confirmed.
2. **Near-equilibrium parts of the state are overshot by a uniform jump.** The thin substellar ice melted back after ×1.5. Scale only what is still spinning up (the tapered ice jump).
3. **Measure profiles and patterns from time means, never from instantaneous restarts.** Two cam.i snapshots gave the right shape for D5's vertical profile but about twice the amplitude (weather).
4. **The slow vertical response of a hot atmosphere is not its fast response.** In D5 the interannual (fast) response warms the whole column by about 1–2 K per K of surface warming. The decades-long (slow) trend cools the mid-troposphere and warms the upper troposphere about twice as fast as the surface. A jump that imposes the fast pattern will relax back.
5. **Judge convergence from energy_top over several years, never from one year or from energy_bot.** Single-year energy_top in D2 ranges from −18 to +6 W/m²; energy_bot is noisier still.
6. **Advice from a post-jump window is contaminated by the post-jump transient.** The ice runs needed about 20 years to settle; a second jump advised from 18 post-jump years was unreliable.
7. **Interim checks need the run directory.** History files reach the archive only when a job ends, so a mid-run check must build the post-jump years from the run directory and join them to the archived series (only the native columns are valid across the join).
8. **Pick probe cases with low internal variability.** D2 has multi-year TS swings of about ±1 K and single-year energy_top excursions near −18 W/m² (year 68: an eight-month shortwave and longwave event, not a data error). A probe response of a few K cannot be read against that in 10 years; D2 needed 20 years and the user's reading of the Gregory plot (see its entry). Its column cloud water (TGCLDLWP, TGCLDIWP) oscillates with a 5–7 yr period, a candidate source of the swings.
9. **In a hot run, energy_top before a jump includes the vapor column still filling.** D4 had energy_top +3 to +6 W/m² and ran about 6.5 W/m² above energy_bot; after a coupled jump that added 17 % vapor, energy_top fell to about 0 and stayed there for 9 years (TS flat at 373.5 K), and TMQ stopped growing. The advisor's local-curve TS_eq (378 K) overstated the equilibrium; a plain annual N-on-TS line over the last 15–20 pre-jump years (373.1–373.6 ± 0.4 K) got it right.
10. **The coupled cam.r jump is numerically safe in CESM 1.2.1 when TEOUT is recomputed.** D4's first-step energy-fixer correction was −1.4e-4 W/m², with no NaNs, and the model ran on without a test segment. S4 and S5, each with +47 % vapor, also completed their first 10 post-jump years without a crash.
11. **Large-imbalance cases are where jumps pay (user, 2026-10-04).** D4 skipped its slow vapor build-up in one coupled jump. Cases with energy_top above about +5 W/m² (S4, S5) are where the real savings are; an advisor refusal for an unresolved slope is not a reason to wait when the slope is small and the imbalance large (a small λ puts the equilibrium further away, not closer).

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

