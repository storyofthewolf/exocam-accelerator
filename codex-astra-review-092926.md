# Code and physics review — 2026-09-29

Reviewed `main` at `eba84bc` (`Add tapered ice jump: scale only conduction-limited ice`). Review only: no implementation changes or live case interventions. The working tree was initially clean; the branch already had one unpublished commit. No project `AGENTS.md` was present; `CLAUDE.md` and the ExoCAM analysis guidance informed the review. In particular, scientific diagnostics must preserve their area weighting and calendar meaning.

The separation between pure numerical code, restart I/O, and case safety is useful. Several earlier issues have been addressed: archive availability is now mandatory by default, missing ICEFRAC/qi columns are rejected, rejected temperature fits block advice, and some input domains are validated. Nevertheless, the current implementation still has failures that undermine its safety claims. I would address findings 1–5 before further production interventions and findings 6–10 before treating the validation metrics as reliable.

## Verification

- Existing suite: **207 passed**, with one viewer test blocked by the sandbox's local socket restriction. That test was rerun with approved socket access and **passed**. Thus all 208 existing tests passed across the two runs.
- Additional synthetic probes reproduced the failures reported below using the existing `tests/synth.py`, `test_check.py`, and restart fixture helpers. Only temporary synthetic NetCDF files were modified; no real restart payloads were read or changed.
- Reviewed the numerical, advisor, taper, checker, restart, rollback, CLI, trend adapter, hindcast, and viewer paths. Cross-checked selected physics claims against local CESM/CICE source. No cluster was contacted, and no coupled-model validation was performed.
- Prior review notes were used to identify previously raised concerns; current code was inspected afresh rather than assuming those concerns remained unfixed.

## Findings

### 1. High — an incomplete archive passes preflight and can produce a false successful rollback

**Locations:** `src/exocam_accelerate/runstate.py:195–209`, `298–312`, `361–366`.

Archive validation checks only whether the named CICE restart exists and lacks the jump attribute. It does not verify the archived component files, the archived rpointers, or their targets. Rollback copies whatever files are present and retires the jump log.

**Reproduced:** an archive containing only `case.cice.r.0101-01-01-00000.nc` passed preflight. After advancing the synthetic run to 0106, its rollback plan also passed. Execution reported success while `rpointer.ice` still named `case.cice.r.0106-01-01-00000.nc`. The continuation would therefore not read the restored state.

**Fix:** construct and validate a complete restart manifest before either operation, including auxiliary files such as `cam.rs` that may not be the first rpointer target. Require archived pointers for the expected components, verify every required target/date/case, and check that archive files cannot alias the writable targets. Verify the restored set before retiring its log.

### 2. High — nonfinite scientific data bypass safety gates, including the final PASS decision

**Locations:** `advise.py:386–404`; `check.py:124–129`, `146–182`, `207–209`; `aqua_ice.py:113–137`; `restart.py:220–231`.

Several gates are comparisons of the form `value > tolerance`; NaN makes these comparisons false. The advisor checks column presence without requiring usable values. The checker counts calendar bins without requiring finite observations in those bins. The restart plugin checks signs but not finiteness.

**Reproduced:**

- Setting all `qi_native`/`qi_int2` values to NaN still produced jump advice with no refusal reasons.
- Setting all ICEFRAC values to NaN did the same.
- Setting all post-jump native N and TS values to NaN returned **PASS**.
- A synthetic restart with positive infinity in `vicen` was written and reported as verified. Equality of two infinities is not evidence of a valid restart.

**Fix:** require finite, physically valid inputs and derived metrics at each boundary. Count only complete valid years, require positive finite ice thickness, validate all serialized fit parameters, reject nonfinite restart values before any write, and make an uncomputed required diagnostic prohibit PASS. Do not silently repair an invalid input restart by independently clamping volume and enthalpy.

### 3. High — the temperature requirement is inconsistent between advice and checking

**Locations:** `advise.py:424–444`; `check.py:159–193`, `207–209`.

The advisor silently skips missing TS, although the standard checker requires an accepted TS reference. Conversely, the checker treats `accepted=True` in the saved reference as sufficient even when the actual TS observations or the linear fit are missing. It adds a reason saying TS was not checked, then returns PASS.

**Reproduced:** removing TS from otherwise valid advisor inputs still allowed a jump. Removing TS observations from an otherwise passing post-jump dataset returned **PASS**.

**Fix:** define one shared required-diagnostics contract. Require an accepted, structurally valid TS fit when advising and finite observed TS when checking. Missing observations should cause an explicit incomplete-data verdict; missing mandatory advice should prevent the jump.

### 4. High — the cold-case advisor accepts the wrong physical regime and can recommend thinning

**Locations:** `phase_space.py:246–253`; `advise.py:344–361`, `412–422`.

The hyperbolic gate tests absolute correlation, so either sign of `b` is accepted. Default targeting halves N toward zero without applying the explicit target's directional check. The factor is clipped to `[1/max_ice_factor, max_ice_factor]`, despite the configuration describing a thickening-only protocol. Positive growth and a physically appropriate deficit/conduction slope are not required.

**Reproduced:** using the existing Stefan fixture with `a=-3, b=+20` yielded accepted advice with factor **0.666667**, `N_now_fit=-2.34750`, `N_after=-2.02123`, and **−43.49 years skipped**. This is outside the intended freezing-dominated acceleration regime.

**Fix:** gate explicitly on the supported cold-case regime, including imbalance convention, slope sign, positive observed growth, and improvement toward the supported target. Require a factor at least one; use the same physical target checks for default and explicit targets, including after range clipping. A rejected regime should be refused rather than converted into a thinning experiment.

### 5. High — `restore` bypasses the safeguards that protect the other restart writes

**Locations:** `cli.py:306–311` (`cmd_restore`); `restart.py:239–252`.

The command immediately overwrites CICE from its backup and retires the log. The statement “before resubmitting only” is documentation, not an enforced condition. It does not check job state, pointer consistency, or whether a subsequent segment has advanced other components. With an explicit old path, it can also retire an old jump's log without rolling back the active coupled state.

**Fix:** enforce that this is an immediate undo at the original segment boundary, with the same scheduler and identity checks as `jump`. If the run has advanced, require whole-set rollback. Revalidate the backup before restoring it.

### 6. Medium — the landing test cannot reliably detect unapplied small jumps

**Locations:** `check.py:51`, `124–133`.

The tolerance is 15% of total expected thickness, rather than a fraction of the proposed increment. Even without natural growth, an unapplied factor F passes when `1 - 1/F <= 0.15`, i.e. for F up to approximately **1.1765**. The subsequent law and growth checks also pass for a healthy untouched control.

**Reproduced:** an entirely unjumped synthetic trajectory checked against a factor-1.1 log returned **PASS**, with landing error −8.51%. Small canary jumps and weak effective tapered jumps are particularly affected.

**Fix:** record the actual pre/post-edit area-weighted restart thickness and verify that the model consumed the intended restart. Evaluate the observed increment against a native-growth control prediction with uncertainty. If a small jump cannot be distinguished from native growth, report that it is unverified rather than asserting it landed.

### 7. Medium — the tapered flux prediction substitutes net thickness change for conduction

**Locations:** `taper.py:97–126`; `restart.py:312–331`; `advise.py:486–494`.

The mask uses `(h_now-h_old)/dt`, and the predictor assumes this net growth scales as `1/f`. That requires more than positive net growth. Restart differences can include basal and surface melt, changing cover, snow-to-ice conversion, and transport where enabled. Old ice area is read and discarded, so even the fixed-area assumption is tested only at the final endpoint.

For a simple local balance `g = C/h - M`, with nonzero thickness-independent melt M, thickening gives `g_after = C/(f*h) - M`, whereas the implemented rule gives `g/f = C/(f*h) - M/f`. They differ even with perfectly measured restarts. A fitted global intercept does not identify the separate local contributions needed to weight a nonuniform jump.

**Numerical stress case:** two equal-area cells with growth `[1, -0.9]` and weights `[1, 0]` give `R(1.5)=-2.3333`. The positive denominator check alone does not guarantee a positive conduction-like ratio or good conditioning. The target solver can stop before that crossing; this example demonstrates the missing domain guarantee, not that every generated advice crosses it.

**Physics cross-check:** the local aqua-FV `ice_therm_vertical.F90:3720–3723` separately computes surface energy and the difference between bottom conductive and ocean heat fluxes. The model itself does not equate net thickness tendency to conduction alone.

**Fix:** establish a quantitative applicability gate from thermodynamic/transport diagnostics and stable cover at both endpoints. Reject strongly cancelling growth/melt budgets or explicitly model their separate contributions. Verify baseline intervals contain no prior artificial jump. Calibrate this prediction on independent coupled runs; the two cases used to select the taper thresholds are not an independent validation set.

### 8. Medium — the hindcast fit window includes data later than its prediction origin

**Locations:** `hindcast.py:157–176`, `266–280`.

For annual samples at half-years and an integer origin, the default `time_tol=0.51` includes the following half-year in the fit. The nearest-origin tie selects the preceding half-year as the baseline. Thus the fit can use a sample a full year later than its baseline. `default_origins` generates precisely these integer origins.

**Reproduced:** with `t=0.5,1.5,...`, `X=t`, origin 10, window 5, and horizon 3, changing only `X(10.5)` by +0.1 changed the prediction from **12.5 to 12.532142857**, although the baseline remained **9.5** and the target observation was unchanged.

**Fix:** resolve the actual origin sample first, fit only through that sample, and measure the forecast horizon from its timestamp. Add a regression test asserting that changing any sample after the origin cannot alter the forecast.

### 9. Medium — gate “false alarms” are classified using the wrong counterfactual

**Locations:** `hindcast.py:287–343`, especially `calibrate_gate`.

For refused proposals, calibration uses persistence error to label `refused_good` and `refused_bad`. Whether persistence succeeds does not establish whether the rejected extrapolation would have succeeded. A nearly stationary oscillation can have good persistence but a dangerous fitted extrapolation; a steadily drifting trajectory can have poor persistence and a good extrapolation. Consequently, `false_alarm_rate` does not measure the quantity advertised by the class and README.

**Fix:** in the offline harness, compute the same clipped extrapolation with only the trust gate bypassed and score that counterfactual. Alternatively rename the existing metric as persistence performance and remove the false-alarm interpretation. This change belongs to offline scoring, not the operational acceptance gate.

### 10. Medium — component alignment validates relative month numbers, not absolute dates

**Locations:** `trend_io.py:85–100`, `129–137`, `140–166`.

Files covering different calendar spans can each contain the same relative month axis `1..N`. `merge_trend_files` accepts those axes as matching, and `case_start_year` chooses the first matching filename. A CAM series for years 1–60 and CICE series for years 61–120 can therefore be paired as if concurrent, corrupting the phase-space fit and reported model year without a shape error.

**Fix:** parse and compare the component filename date spans, validate their lengths against the month index, and align on absolute calendar months. Reject ambiguous multiple spans or select a single complete shared span explicitly.

## Additional physics and engineering improvements

- **Albedo preservation is overstated.** `aqua_ice.py:11–16` asserts that unchanged area and surface temperature imply unchanged albedo. Local CESM `ice_shortwave.F90:716–723` computes bare-ice albedo from `vicen/aicen`; radiation transmission also depends on thickness. The statement can be approximately valid for thick ice where albedo has saturated, but it is not guaranteed for all supported factors/cells or optional snow scaling. Preserve the narrow thick-ice rationale, verify the active radiation scheme and affected categories, and assess first-coupling flux changes before claiming coupled consistency. Modern CICE documentation likewise describes albedo being computed consistently with both area and thickness: [CICE fundamental variables](https://cice-consortium-cice.readthedocs.io/en/cice6.0.1/science_guide/sg_fundvars.html).
- **The tapered law offset is held constant indefinitely.** `check.py:148` adds the jump-time offset at every subsequent thickness. Under local Stefan evolution, differently scaled cells evolve at different rates, so a fixed additive correction is an approximation. Bound its validity to a calibrated watch interval or derive a time-dependent spatial reference. This is a scientific validation gap, not a demonstrated error for the existing short canaries.
- **Restart writes are not transactional.** `restart.py:220–234` overwrites variables in place before verifying and writing the log. An interruption can leave a partially edited active restart with no active log for the normal rollback selector. Stage and verify a complete output, maintain a recovery journal, then atomically replace the active file on the same filesystem. A pristine backup helps recovery but does not make a partial write safe to resume.
- **Bind the exact target and case throughout.** Preflight compares pointer dates, not the resolved ice pointer target against an explicit `--cice-r`; same-date files from different cases can evade that check. `cmd_check` also does not compare the requested trend case with the jump log's case. Require exact identities before mutation or scoring.
- **Small prior jumps evade automatic history separation.** `advise.detect_jumps` requires an annual thickness increase above 8%. A tapered jump with a smaller effective increment, or a small canary, can remain in a later fitting window. Use authoritative jump records when available; threshold detection should be a fallback.
- **Repeated I/O and numerical work can be reduced.** `cmd_jump` calls the full restart calculation for preview and again for application. Viewer `/api/cases` computes complete plot payloads for every case, and selecting a case recomputes them. `fit_saturating` performs hundreds of small least-squares fits for each reference. Cache annualized columns and fits by input fingerprint plus configuration; provide a summary-only viewer path. Reuse prepared restart work only while verifying source identity and freshness before committing the write.
- **Reproducibility metadata is incomplete.** `Advice.to_dict` omits settings such as enthalpy tolerance, gate thresholds, and settling/detection controls. Store the full effective configuration and code version with the existing input hashes.
- **Documentation contradicts the implementation.** README still says plugins and restart writing are unimplemented and cites 152 tests, while later sections describe the implemented workflow and the suite contains 208 tests. Update the overview, status, and stale module descriptions together.

## Suggested order

1. Repair archive completeness, restore safeguards, finite-data handling, TS requirements, and cold-regime enforcement.
2. Add focused regression tests for the reproduced failures, exact case/date binding, and small-jump landing detection.
3. Repair hindcast temporal leakage and counterfactual calibration before using those metrics to choose thresholds.
4. Validate taper energetics, evolving spatial structure, and first-coupling behavior on independent canaries and untouched controls.
5. Improve caching, provenance, and documentation after the safety and scientific contracts are explicit.

The existing tests demonstrate that the implemented algebra and intended happy paths work. They do not yet establish that missing diagnostics, incomplete restart sets, or unsupported physical regimes are refused, nor that a global-mean PASS establishes spatial or coupled-model fidelity.
