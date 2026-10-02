# Session handoff, 2026-10-02

Read CLAUDE.md first, then README.md (including the Glossary) and docs/ocean-jump.md. Load hpc-connect before any Discover contact and exocam-rules before any case action. This file replaces all earlier handoffs.

## State at end of session

Both repos are committed on main locally:

- **exocam-accelerate**: f1a9ee6 (atm-profile from time means) plus the commit that adds this file.
- **exocam-trend** (`../exocam-trend`): 0b0ee7e (per-level global means, `--profile`).

**Neither is pushed yet.** The agent's `git push` was refused by the permission classifier, so the user pushes both (`git push origin main` in each). The Discover checkouts (`~/pythonWorkspace/exocam-accelerate` at 644f6e9 and `~/pythonWorkspace/exocam-trend`) are behind until the user has pushed and they have been `git pull`ed through hpc-connect. On Discover the package is not pip-installed: run it as `PYTHONPATH=~/pythonWorkspace/exocam-accelerate/src python -m exocam_accelerate ...` with `export PATH="$HPC_PYTHON_BIN:$PATH"`.

## What changed this session

- **exocam-trend 0b0ee7e:** `trend.py --profile T,Q` reduces 3D cam.h0 fields to area-weighted global means on each **model level** (no interpolation; user decision), in the same pass over the files as the 2D fields. With `--save-data` it writes `data/<case>_<first>-<last>_camlev_<VAR>.txt` per field plus `_camlev_PMID.txt` (hyam*P0 + hybm*<PS>), header `month L01 ... Lnn` (L01 = top), native monthly means only. `run_trend_batch.sh --profile T,Q` passes it through. `--outdir` now copies only the files that run wrote; it used to copy every older series of the case, and `load_case` then merged them (that broke the D5 advisor this morning). README/CLAUDE.md there now name the real `--save-data` flag.
- **exocam-accelerate f1a9ee6:** `atm-profile` measures the coupled jump's vertical profile from **time means**, never from cam.i (user decision: instantaneous restarts carry weather). It reads `--trend-dir` (the `_cam.txt` for TS plus `_camlev_T.txt` and `_camlev_PMID.txt`), builds annual means, and with the window ending at the advice's model year uses `--method trend` (default: ratio of least-squares slopes of each level and TS over `--window 10` years, delta-method standard error per level) or `--method annual` (difference of two annual means). `--archive`, `--rundir` and `--domain-file` are gone from atm-profile, and the cam.i code (`build_atm_profile`, `_cam_i_profile`, `archived_cam_i`) is removed. New: `atmos.profile_from_annual`, `trend_io.level_annual_means`; the profile file and viewer carry `raw_gain_se`. 297 tests pass. The jump still edits the instantaneous cam.r/docn.r as before.

## Why: D5 year-60 findings (2026-10-02)

D5 finished at restart 0061 (archive and rpointers checked). Trends through year 60 are in `$HPC_SCRATCH/atlasfu_d5_jump/trend60/` (only the 0060-12 `_cam.txt`). `advise-ocean` (defaults: local curve, energy_top, 40-yr window) gave TS 376.07 K, energy_top +2.75 ± 0.31, alpha_diff 2.01 ± 0.23 (year 50: 1.12 ± 0.11), TS_eq 377.7 K (year 50: 380.8), recommended TS step +0.83 K (year 50: +2.9), ocean heat fraction 0.27 (`$HPC_SCRATCH/atlasfu_d5_jump/D5_advice60.json`).

The old cam.i atm-profile (0051 to 0061) gave -6 K/K through 400-1400 hPa and +20 K/K at the model top, so the auto top cut the jump off at 3927 hPa (mass-weighted gain 0.09). The user identified this as weather in instantaneous fields; that led to the redesign. `D5_coupled60.json` and its `.atmprofile.nc` in that directory are from the old method: **do not use them for a jump.** Local copies are in `../scratch/atlasfu/d5view/`.

Cost measured on Discover (reading only T, Q and PS; cam.h0 is uncompressed netCDF-3, 45.7 MB per month, 0.68 MB per 3D field): two annual means 3 s, a 10-yr trend 11 s. File opens dominate. The "~5 GB per 10 yr" only matters if whole files are copied or processed.

## User decisions to respect

- energy_top (TOA) is the most important convergence quantity. Never judge convergence from energy_bot or from a single year.
- The Gregory relation need not be linear for hot planets.
- Curved or noisy trajectories are sometimes worth jumping anyway; the override exists for that.
- Runaway is a warning; the rollback call is the user's.
- The atmosphere jump never touches p < 100 hPa.
- **The profile comes from time means of a prior window, never from instantaneous restarts. Model levels, no interpolation.** Post hoc the user wants to discuss different window lengths.
- **Horizontal layer means stay** (lon-lat gradients are small in the hot lower atmosphere). Revisit 2D spatial jumps only if spin-up rates differ by location (e.g. sea ice), as the ice taper does.
- The D5 jump is IN PLACE (no clone, no control run).
- No automatic polling of cluster jobs; the user says when a job has finished.
- Work on main in both repos (user, 2026-10-02).

## Runs in flight

**D2 (`exocam_atlasfu_D2`).** Ocean-only probe at restart 0061 (somtp +8 K uniform, f_ocean 0.5; advice `$HPC_SCRATCH/atlasfu_d2_probe/D2_probe4_hr2.json`). Continuation job 58680001 (years 61-70, RESUBMIT=0) was RUNNING at 13 h elapsed on 2026-10-02 morning. When the user says it has finished: regenerate trends into a fresh directory (`run_trend_batch.sh --cam --nmonths 840 --int1 1 --int2 10 --outdir <fresh dir> exocam_atlasfu_D2`), then `check`. Expect WAIT and a 5-10 yr extension.

**D5 (`exocam_atlasfu_D5`).** Finished at 0061, not queued, untouched. Nothing has been jumped.

## NEXT TASK: D5 profile from time means, then decide the jump

1. After the user has pushed: `git pull` both repos on Discover (hpc-connect; check `bin/hpc-sync-status`).
2. Regenerate D5 trends with profiles into a fresh directory: `cd ~/pythonWorkspace/exocam-trend && ./run_trend_batch.sh --cam --profile T,Q --nmonths 720 --int1 1 --int2 10 --python python --outdir $HPC_SCRATCH/atlasfu_d5_jump/trend60p exocam_atlasfu_D5`. That needs `--allow-exec --allow-write` with the user's OK; about 1-2 min. The directory should then hold exactly one `_cam.txt` and three `_camlev_*` files.
3. `advise-ocean` on that directory (it ignores the camlev files) and confirm the year-60 numbers above. Then run `atm-profile --advice ... --trend-dir ... --json ...` with `--method trend` and `--method annual`, and with a few `--window` values (5, 10, 20) for the post-hoc discussion. These write to scratch only.
4. Pull the outputs into `../scratch/atlasfu/d5view/` and relaunch the viewer (`exocam-accelerate view ../scratch/atlasfu/d5view`). It runs as a background command that hits a time limit, so start it with a long timeout. Report the profiles, their standard errors and the mass-weighted gain to the user.
5. Only with the user's decision: run the earlier launch sequence. `jump --dry-run` (check the vapor increase is well under the 50 % refusal), then `jump` in place (pristine backups of cam.r and docn.r go to `run/exocam_accelerate/`), then ONE 1-month segment (`runmgr.py continue exocam_atlasfu_D5 --set STOP_OPTION=nmonths --set STOP_N=1 --set RESUBMIT=0`, preview then `--execute`). Inspect `atm.log` when the user says it is done: energy-fixer first-step correction vs D5's normal restarts, no NaNs, TS and PS change. STOP and report before any 10-yr continuation. Roll back with `exocam-accelerate rollback` if anything is wrong.

The coupled atmosphere jump has been rehearsed on D4 copies but never run in-model.

## Open questions and later

- Whether to extend D2 after its first check.
- The D4 profile numbers in docs/ocean-jump.md come from the old cam.i method; re-measure if D4 is jumped.
- Optional diagnostic: zonal means per level, to show when horizontal layer means stop being adequate.
- The user may want TS_now as a 5-yr mean too (agent suggestion, not decided).
- The viewer's profile ±se band was not browser-tested.
