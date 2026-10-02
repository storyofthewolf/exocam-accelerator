# Session handoff, 2026-10-01

Read CLAUDE.md first, then README.md (including the Glossary) and docs/ocean-jump.md. Load hpc-connect before any Discover contact and exocam-rules before any case action. This file replaces all earlier handoffs.

## State at end of session

Local, GitHub (origin main) and the Discover checkout (`~/pythonWorkspace/exocam-accelerate`) are all at the same commit: the one that adds this file (`git log -1 -- docs/session-handoff.md`; its parent is 694ddce). There are no feature branches.

Push with `git push origin main` (works from Claude). Discover updates by `git pull` through the hpc-connect skill. On Discover the package is not pip-installed: run it as `PYTHONPATH=~/pythonWorkspace/exocam-accelerate/src python -m exocam_accelerate ...`.

## What changed this session

- a5cd2c9 (with 9f34570): feature/ocean-jump merged to main; the feature branches are deleted.
- ccf3872: advisor `N_now` is the 5-yr mean with its standard error; the probe heat fraction uses the longest window.
- 9d45322: probe check gains a TS-trajectory readout (tau, TS_eq, side, lambda = C/tau) alongside energy_bot.
- dcd2df5: probe check reads energy on energy_top (energy_bot is information only), noise from a long detrended pre-window, 3 sigma; energy false-alarm FAILs downgraded.
- 34ddfbd: "runaway greenhouse suspected" is a WARNING with verdict WAIT, never an automatic FAIL (user decision).
- 3d45ce5: atmosphere jump is troposphere-only: nothing at p < 100 hPa, log-p taper from 200 hPa (`--ceiling`, `--taper-bottom`).
- 1843270: "heat ratio" is presented as the ocean heat fraction f_ocean = C_ocean/C_total (`--ocean-fraction`; `--heat-ratio` kept as a hidden old name). Also fixed the viewer slider bug (it always sent heat_ratio=1).
- 3efc25c: advisor default is the local curve: a quadratic N(TS) about the current state on native annual means, alpha_diff = -dN/dTS (Gregory 2004; Wolf et al. 2018 JGR-A), on energy_top. The time-domain gate passes decelerating curvature and needs 2 sigma for a sign flip. `--override-gate` added.
- b699abb, 58a0be3, 694ddce: advise-ocean printout names the local fit's series; README glossary and ocean section updated.

## User decisions to respect

- energy_top (TOA) is the most important convergence quantity. Never judge convergence from energy_bot or from a single year.
- The Gregory relation need not be linear for hot planets.
- Curved or noisy trajectories are sometimes worth jumping anyway; the override exists for that.
- Runaway is a warning; the rollback call is the user's.
- The atmosphere jump never touches p < 100 hPa.
- The D5 jump is IN PLACE (no clone, no control run).
- No automatic polling of cluster jobs; the user says when a job has finished.

## Runs in flight

**D2 (`exocam_atlasfu_D2`).** Ocean-only probe applied in place at restart 0061: somtp +8.000 K uniform (area mean 351.885 to 359.885 K), explicit `--probe-dt 4 --heat-ratio 2` (f_ocean 0.5). Advice: `$HPC_SCRATCH/atlasfu_d2_probe/D2_probe4_hr2.json`; backup and `.accel.json` are in `run/exocam_accelerate/`. Continuation job 58680001 (STOP_N=10 nyears, RESUBMIT=0) covers model years 61-70.

D2 context: it was NOT converged at year 60 (20-yr Etop 3.5 +/- 1.2 W/m2, TS +0.16 +/- 0.04 K/yr), and the local-curve alpha_diff was 0.65 +/- 1.15 (unresolved), which is why a probe was used.

Next, when the job completes: regenerate trends (exocam-trend `run_trend_batch.sh --cam --nmonths 840 --int1 1 --int2 10 --outdir $HPC_SCRATCH/atlasfu_trend exocam_atlasfu_D2`), then run `check`. Expect WAIT (synthetic runs gave PASS in about 43 % of cases at 8 settled years and 96 % at 13), so extend by 5-10 yr.

**D5 (`exocam_atlasfu_D5`).** The last segment job 58668877 (RESUBMIT=0) ends at restart set 0061, expected early 2026-10-02. The archived 0051 set has cam.i, cam.r and docn.r.

## NEXT TASK: D5 coupled-jump launch sequence

User-approved. To be executed by a Sonnet agent with a scoped D5 grant once the user confirms D5 has finished.

1. Confirm D5 is not queued or running (`bin/hpc-jobs --name exocam_atlasfu_D5`). Confirm the archive has `rest/0061-01-01-00000` with cam.i, cam.r and docn.r, and that the rpointers name 0061.
2. Regenerate D5 trends through year 60 (`--nmonths 720`) into a fresh directory containing only that file.
3. Run `advise-ocean` with defaults (local curve, energy_top), full advisor step (n_fraction 0.5), `--json`. At year 50 it advised TS 374.9 to 377.8 K (alpha_diff 1.12 +/- 0.11, TS_eq about 380.8 K). Recompute on year-60 data; stop and report if it refuses or differs wildly.
4. Run `atm-profile --advice ... --archive $HPC_ARCHIVE/exocam_atlasfu_D5 --rundir RUN --json coupled.json` (profile from cam.i 0051 to 0061; troposphere-only defaults). In the coupled jump somtp moves by delta-TS itself (ocean heat fraction 1), atmosphere T by gain times delta-TS, q at fixed RH, dry mass kept, TEOUT recomputed. Constants come from the run's atm.log, and the tool refuses unless TEOUT is reproduced to 0.1 %. Check that the vapor increase is well under the 50 % refusal (expect about 10 %).
5. Run `jump --dry-run`, then `jump` (in place; pristine backups of cam.r and docn.r go to `run/exocam_accelerate/`). Verify the edits.
6. Submit ONE 1-month segment: `runmgr.py continue exocam_atlasfu_D5 --set STOP_OPTION=nmonths --set STOP_N=1 --set RESUBMIT=0` (preview first, then `--execute`).
7. When that month completes (the user will say; no polling), inspect `atm.log`: the energy fixer's first-step correction compared with the same lines from D5's normal restarts (should be its usual size), no NaNs or crash, the month's global TS (about +2.9 K vs pre-jump) and PS (a rise of order tens of hPa).
8. STOP and report to the user before the 10-year continuation (STOP_OPTION=nyears STOP_N=10). If anything is wrong, run `exocam-accelerate rollback` to the pristine 0061 set (restores both files).

The coupled atmosphere jump has been rehearsed on D4 copies but never run in-model; this is its first real test.

## Open questions and later

- Whether to extend D2 after its first check.
- The user may want TS_now as a 5-yr mean too (agent suggestion, not decided).
- Viewer changes this session were not browser-tested.
- The 3-bar and D1-D3 hindcast numbers are in docs/ocean-jump.md.
