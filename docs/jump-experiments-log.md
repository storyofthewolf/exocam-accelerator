# Jump experiments: what worked, what did not

A running log of every in-model jump and every measurement that changed the method. Add a dated entry for each experiment (newest at the bottom of its section), and update **Lessons so far** when an entry changes what we would do next time. Record negative results as carefully as positive ones. Design detail lives in `docs/ocean-jump.md` and `docs/restart-integration-questions.md`; this file records outcomes.

Entry template: case, restart, what was changed and by how much, the advice it came from, what we expected, what happened (numbers), verdict, lesson.

## Lessons so far

1. **A jump keeps only the part that lies along the slow trajectory.** Sea ice (a slow variable on its own) held almost all of a ×1.5 jump. An ocean-only somtp jump in a hot run (D2) lost most of its step within 3 years even though the whole atmospheric column followed the ocean within the first year: the perturbation is pulled back by the fast feedback (λ ≈ 2.4 W/m²/K), while the run's own slow drift continues at its old rate.
2. **Near-equilibrium parts of the state are overshot by a uniform jump.** The thin substellar ice melted back after ×1.5. Scale only what is still spinning up (the tapered ice jump).
3. **Measure profiles and patterns from time means, never from instantaneous restarts.** Two cam.i snapshots gave the right shape for D5's vertical profile but about twice the amplitude (weather).
4. **The slow vertical response of a hot atmosphere is not its fast response.** In D5 the interannual (fast) response warms the whole column by about 1–2 K per K of surface warming. The decades-long (slow) trend cools the mid-troposphere and warms the upper troposphere about twice as fast as the surface. A jump that imposes the fast pattern will relax back.
5. **Judge convergence from energy_top over several years, never from one year or from energy_bot.** Single-year energy_top in D2 ranges from −18 to +6 W/m²; energy_bot is noisier still.
6. **Advice from a post-jump window is contaminated by the post-jump transient.** The ice runs needed about 20 years to settle; a second jump advised from 18 post-jump years was unreliable.
7. **Interim checks need the run directory.** History files reach the archive only when a job ends, so a mid-run check must build the post-jump years from the run directory and join them to the archived series (only the native columns are valid across the join).
8. **Pick probe cases with low internal variability.** D2 has multi-year TS swings of about ±1 K and single-year energy_top excursions near −18 W/m² (year 68: an eight-month shortwave and longwave event, not a data error). A probe response of a few K cannot be read against that in 10 years.

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

## Hot, ice-free ocean (atlasfu D-series, 4 bar, 340–375 K)

### 2026-10-01: D2, ocean-only probe, docn.r somtp +8 K at 0061 — WAIT, little gained

- **Change:** somtp +8.000 K uniform in place (351.885 → 359.885 K), probe mode, probe-dt 4, assumed ocean heat fraction 0.5, so a TS step of +4 K. Advice `$HPC_SCRATCH/atlasfu_d2_probe/D2_probe4_hr2.json`. Continuation job 58680001, years 61–70; archived through restart 0071.
- **Final check (2026-10-02, years 61–70, 8 settled):** verdict WAIT.
  - Landed: first-year TS +5.33 K net of drift against +4.00 expected, so the implied ocean heat fraction was 0.67, not 0.5.
  - TS: 357.87 (61), 356.78, 354.92, 353.79, 354.01, 354.24, 355.42, 355.61, 355.38, 353.80 K (70). Pre-jump drift +0.163 ± 0.054 K/yr; post-jump (years 63–70) +0.07 ± 0.22.
  - Against the pre-jump trend extrapolated (years 41–60 fit), D2 averages +0.8 K ahead over years 64–70 (about 5 years of spin-up), but year 70 is 0.5 K *behind*. With multi-year swings of about ±1 K, the net gain is between zero and about 1 K: not distinguishable from zero with confidence.
  - energy_top +3.4 → −2.2 ± 3.3 W/m², λ = 2.4 ± 1.4 W/m²/K: restoring, no sign of runaway, not resolved. Year 68 alone has energy_top −18.1 W/m²; without it the post-jump mean is +0.05.
  - Year 68 is a real event, not a data error: energy_top −5 to −30 W/m² in every month, first from higher OLR (FLNT +15 to +31 W/m² in Mar–Jul), then from less absorbed sunlight (FSNT −16 to −44 W/m² in Aug–Nov), with TS about +1 K throughout. TS then fell 1.6 K in year 70. D2 had a similar 1.5 K dip in year 58 before the jump.
  - The atmosphere followed the ocean: in year 61 the whole column warmed (+15 K at 144 hPa, +3 K at 660 hPa, +10 K at 2850 hPa; 0.5–2.7 K per K of TS), and over years 64–70 it stays 3–5 K warmer than years 56–60 at every level, with water vapor up 10–95 %. D2 (TS ≈ 355 K) is in the uniform-warming regime that D5 left near TS ≈ 372 K.
- **Interpretation (not proven):** the jump is restored by the fast feedback (λ ≈ 2.4 W/m²/K, relaxation within about 3 years), while the run's slow drift, with N ≈ +3.5 W/m² sustained for decades and a flat Gregory relation, continues at its own rate. The fast λ restoring a perturbation is not the slow effective λ along the trajectory. A state jump only helps if it moves whatever sets the slow drift, which we have not yet identified for D2.
- **Lesson:** the ocean-only probe was safe (no crash, no runaway) but bought little. D2 is a poor probe target because of its internal variability.
- **Decision (2026-10-02, user):** extend D2 by 10 years in place (job 58691769, years 71–80) to decide whether ocean-only somtp jumps are worth doing at all. Rollback not indicated.

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
