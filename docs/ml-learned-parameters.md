# Extension: learning the relaxation parameters across a suite of runs

Design note (2026-07-17), companion to `tier0-hindcast-findings.md` §8. **Nothing
here is implemented or committed as code.** This records a direction and, more
urgently, one thing worth doing *now* regardless of whether the ML is ever built.

## The idea

> **Update (2026-07-18):** `phase-space-extrapolation.md` supersedes the
> time-domain framing below — the state is better extrapolated against the energy
> imbalance `N` to `N = 0` than against time to `t → ∞`. The learned target
> therefore shifts from "infer `τ`" to "infer the phase-space feedback relation
> `Ts(N)` (its slope `λ` and any curvature/kink) from a short burst + case
> descriptors." The data roadmap and small-data/GP reasoning below are unchanged
> and apply verbatim to the phase-space parameters.

Tier-0 §8 established the step model: the approach to equilibrium is a saturating
relaxation, and stepping along a fitted

    X(t) = X_eq − A · exp(−(t − t0) / τ)

beats a linear (or naive-quadratic) Euler step on every variable at every Δt. The
weakness of that step is not the *form* — it is the *fit*. Each `(X_eq, A, τ)` is
estimated from one short trailing window, which is noisy: a 30-yr window early in
a 200-yr relaxation barely samples the curvature, so `τ` and the asymptote are
poorly constrained. This is exactly why the single exponential under-predicts
`hi` at large Δt (§8, Fig 11).

**The proposal:** rather than fit `(X_eq, A, τ)` fresh from each window in
isolation, *learn* a model that infers them from (a) features of the observed
early trajectory and (b) static case descriptors — instellation, obliquity,
rotation period, initial ice fraction, whatever distinguishes the runs — trained
across a *suite* of past relaxations. The physics fixes the curve family; the
learned model supplies better parameters than a local fit can, and can predict
the equilibration timescale of a *new* run from its first few decades.

## Why this is well-posed here (and where it is not)

**This is a small-data, strong-prior problem, not a big-model problem.** The
binding constraint is training signal, not model capacity. The number of *runs*
— not the number of monthly samples — sets the effective sample size, because
windows drawn from one run are heavily autocorrelated (near-duplicates). A
flexible neural network will overfit a suite of tens of runs. So:

- **Prefer few interpretable parameters over a free function approximator.** The
  target is 2–3 numbers per variable (`X_eq`, `τ`, and `A` which is then pinned
  by the observed state), not a trajectory. Predict *those*.
- **Prefer Gaussian-process regression to an NN at this scale.** A GP over the
  relaxation parameters gives calibrated predictive *uncertainty* for free — and
  the trustworthiness gate can consume that uncertainty directly (refuse a step
  when the inferred `τ`/`X_eq` is too uncertain). An NN gives a point estimate the
  gate cannot interrogate. Move to an NN only when the suite is large enough
  (many hundreds of runs) that the GP's cost or rigidity bites — see the data
  roadmap below.
- **The coupled-mode structure is a gift.** Tier-0 §7 showed TS, TOA energy
  balance, and the ice reservoir move as one eigenmode (|r| > 0.9). The per-run
  dynamics therefore live on a low-dimensional manifold, so a model conditioned
  on a few case descriptors has a real chance of generalizing `τ` to unseen runs
  — the whole point.

**What to avoid.** A model that ingests raw time series and emits the accelerated
restart state *directly* discards the conservation structure (energy balance, the
coupled mode) and cannot be audited — fatal for a tool whose entire premise is a
refuse-when-untrustworthy gate. A black-box jump across a nonlinear feedback
threshold is precisely the failure the gate exists to prevent. Any ML component
must (1) predict *parameters of the physical curve*, not the state, and (2) emit
an uncertainty the gate can threshold. Those two rules keep the extension inside
the existing safety architecture instead of bypassing it.

## Tiers (each strictly beats the one below; stop when the win plateaus)

- **Tier A — analytic fit, no ML (already in §8).** `curve_fit` the saturating
  form per run/window. This is the baseline any ML must beat — not "beat linear."
- **Tier B — learn a prior over the fit parameters (the sweet spot).** A GP maps
  (early-trajectory features + static case descriptors) → `(X_eq, τ)` with
  uncertainty. Replaces the noisy per-window fit with a suite-informed estimate;
  predicts a new run's timescale from its opening decades. This is the direct
  realization of "infer the appropriate curvature parameters."
- **Tier C — learn the residual / a second timescale.** The single exponential
  under-predicts `hi` at large Δt, likely because ice relaxes on a slower
  timescale than the atmosphere. Let a small model predict the *correction* to
  the analytic fit (or a two-timescale `τ1, τ2`). Because it predicts a residual
  around a physical fit, it cannot run away the way the raw quadratic did.

## The data roadmap — the actually strategic part

You have three data sources with very different economics, and the right move
depends on treating them differently:

1. **Legacy time series (available now).** The seed training set. First job:
   inventory them and normalize to the exocam-trend text schema this tooling
   already parses (`trend_io`), recording the static descriptors (instellation,
   obliquity, rotation, initial state) as per-run metadata — those descriptors
   are the model's input features and are useless if not captured alongside the
   series.
2. **Purpose-generated runs (available but expensive).** Do **not** spend these
   on random coverage. Once a Tier-B model exists, use *active learning*: fit the
   GP, find where its parameter uncertainty is largest (which corner of
   instellation/obliquity/ice space it cannot yet predict), and generate a run
   *there*. Each expensive run then buys the most model improvement possible.
   Until the model exists, prioritize spread across the descriptor axes you most
   expect to matter (instellation, initial ice), not repeats near existing runs.
3. **Every future production spin-up (free).** This is the compounding source and
   the reason to act now. **Any** ExoCAM run to equilibrium is a complete
   relaxation trajectory — free training data — *if the time series is captured
   as it runs*. The single highest-leverage action available today, well before
   any ML is built, is to make trajectory capture a standing side-effect of
   normal runs: emit the exocam-trend series (plus the run's static descriptors)
   automatically for every spin-up and archive it to a common store. The training
   set then grows for free with normal science, and by the time a Tier-B model is
   worth training the data is already there instead of being reconstructed after
   the fact.

**Recommendation, independent of the ML timeline:** stand up the trajectory
archive and the automatic capture now. It is cheap, it is decoupled from every
open design question in the main findings, and it is the prerequisite that has a
lead time measured in runs. The model can wait; the data collection should not.

## Validation discipline (stated up front so it is not forgotten)

- **Leave-one-run-out, never a random split.** Windows from the same run are
  near-duplicates; a random train/test split leaks them across the boundary and
  flatters the model badly. The only honest test of "can it predict a new run" is
  to hold out entire runs.
- **The metric is downstream hindcast error, not parameter MSE.** Score a learned
  `(X_eq, τ)` by the Tier-0 hindcast error of the step it produces (the §8
  harness already computes this), not by how close the parameters are to a local
  fit — a parameter set can be "wrong" yet yield a better step, and vice versa.
- **Regression guard against the analytic baseline.** The learned model ships only
  if it beats the Tier-A `curve_fit` step under leave-one-run-out. If it does not,
  the analytic fit stands and the ML is not worth its complexity.

## Where this sits relative to the main findings

This is downstream of, and gated by, `tier0-hindcast-findings.md` §9 item 1: the
saturating-relaxation step must exist as the core stepper before there is a
parameter for ML to improve. The sensible order is (1) implement the analytic
saturating step + its guard (Tier A), (2) stand up trajectory capture and the
archive (independent, do early), (3) revisit Tier B once the suite has grown
enough to train and validate leave-one-run-out.
