# Phase-space extrapolation: energy imbalance as the governing coordinate

Session note (2026-07-18), offline bake-off. Companion to
`tier0-hindcast-findings.md` §8 (the time-domain saturating step) and
`ml-learned-parameters.md`. **Offline analysis only; no core code changed.**
Driver: `../scratch/phase_space_bakeoff.py`.

## The idea being tested

Tier-0 §8 extrapolated the state against **time** (`Ts(t) = X_eq − A·e^(−t/τ)`,
step to `t → ∞`). That requires estimating `τ` and a time-asymptote from a short
window — ill-conditioned early (τ barely sampled) and, as shown below, unstable
near equilibrium (a flat window cannot constrain τ at all).

The alternative, motivated by the user: drop time as the independent variable and
work in **phase space**, with the net energy imbalance `N` (TOA, here
`energy_top`, W/m²) as the governing coordinate. Equilibrium is the physical
condition **`N = 0`**. Instead of extrapolating to `t → ∞`, fit the state as a
function of `N` and evaluate at `N = 0` — a *known x-value*, pinned by physics.

This is inspired by Gregory-style ΔN–ΔT regression (inferring equilibrium from a
short transient), but with the deliberate generalization that **the fit need not
be linear** — a saturating or otherwise curved `Ts(N)` may predict the endpoint
better when feedbacks are state-dependent (ice-albedo, deep-atmosphere
adjustment).

## The bake-off

Three predictors, each fitting only a window and predicting `Ts_eq` (= actual
final annual-mean Ts, the best available equilibrium proxy; these runs end with
residual `N ≈ −1.3` W/m²):

| predictor | form | endpoint |
|-----------|------|----------|
| `time-exp` | `Ts(t) = X_eq − A·e^(−t/τ)` (time domain, §8) | `X_eq` |
| `gregory`  | linear `Ts = c0 + c1·N` | `Ts(N=0) = c0` |
| `phase-nl` | curved `Ts(N) = T_eq − A·(1 − e^(N/N0))` | `Ts(0) = T_eq` |

Three fit-window regimes: **near-eq** (late 50 yr — the *clean end-state* jump),
**early** (yr 10–50 — already near N=0), **wide** (yr 2–18 — spans the
fast, high-imbalance early branch; the "jump from far out" test). Median absolute
`Ts_eq` error over the 15 cases, `int2` fit:

| regime | time-exp | gregory | phase-nl |
|--------|---------:|--------:|---------:|
| near-eq | 0.42 (mean **19.9**) | 0.92 | 1.08 |
| early   | 1.12 | **0.69** | 1.00 |
| wide    | 1.74 (mean **1.5e3**) | 0.29 | **0.15** |

(means in **bold** where they blow up: time-exp `max` reached 207 K near-eq,
2.3e4 K wide.)

## Findings

**1. Phase space beats the time domain decisively — the endpoint being pinned is
the whole point.** `time-exp` is not merely worse but *unstable*: when the fit
window is flat (near equilibrium) or short (wide), `τ` is unconstrained and `X_eq`
flies off (mean error 20 K near-eq, 1500 K wide, max 2.3e4 K). `gregory` and
`phase-nl` stay bounded and small in **every** regime, because extrapolating to
the fixed `N = 0` cannot run away the way extrapolating to `t → ∞` can. This is a
quantitative vindication of the premise: **energy imbalance is the more
fundamental coordinate than time** — the equilibrium condition lives on the
N-axis, not at infinite time. It also means the *clean end-state* use case (jump
from a nearly-converged run to true `N=0`) should use the phase-space fit, not the
time-domain one.

**2. Nonlinear beats linear where the physics is nonlinear — but only on smoothed
input.** On the wide window (the widest N-range, where curvature is real and
resolved), `phase-nl` (0.15 K) beats `gregory` (0.29 K) ~2×. Geometrically
(Fig 12): for most cases `Ts(N)` is strikingly linear and both fits nail the
`N=0` star; but **pt10 has a visibly concave `Ts(N)`** — there the linear Gregory
overshoots the endpoint while the curved fit bends with the true relation and
lands close. That is the ice-albedo feedback's state-dependence showing up as
phase-space curvature, exactly the case the linear special-case cannot handle.
This is the empirical support for *not limiting the method to linear regression*.

**3. The nonlinear advantage is conditional on smoothing.** With `native`
(unsmoothed) input, `gregory` (0.61) beats `phase-nl` (0.84) on the wide window —
the noisy fast early drop lets the curved form overfit spurious curvature. Only
with `int2` (long smoothing) does clean curvature emerge and `phase-nl` win. This
ties the phase-space direction to the longer-smoothing hypothesis
(`tier0-hindcast-findings.md` §4): **long smoothing is a prerequisite for
nonlinear phase-space extrapolation**, not just a nicety.

**4. In these cold cases the huge imbalance collapses fast** (N: −40 W/m² at
yr 5 → −5 by yr 20), so "jump from year 10" is, in phase-space terms, a *short*
extrapolation already near N=0 — which is why all methods do tolerably on the
early window. The genuinely long extrapolation is the wide (yr 2–18) window that
reaches back to the high-N branch; that is the real test of the asynchronous
ladder's first jump, and phase-space handles it (median 0.15–0.29 K) where time
domain fails.

![phase-space bake-off](figures/tier0/fig12_phase_space_bakeoff.png)

*Fig 12 — `Ts` vs net imbalance `N`, `int2` fit. Blue = wide fit window (yr 2–18);
red dashed = linear Gregory; green = nonlinear `Ts(N)`; star = true equilibrium at
`N=0`. pt01/pt13 are near-linear (both fits succeed); pt10 is concave — Gregory
overshoots, the curved fit tracks.*

## Consequences for the method

- **Adopt phase-space (`state`-vs-`N`) as the extrapolation frame**, with `N = 0`
  the endpoint. The time-domain saturating step (§8) is its shadow — same physics,
  worse conditioning — and should be demoted to a fallback where an imbalance
  series is unavailable.
- **The fit family is pluggable and problem-class-specific.** Linear Gregory is
  the near-equilibrium / constant-feedback special case; a curved `Ts(N)` wins
  under state-dependent feedback (cold ice worlds); other classes (Earth-twins,
  mini-Neptunes) will want their own forms. This mirrors the existing
  `VariableConstraint` plugin pattern, one level up: **register the phase-space
  fit form per problem class.**
- **The gate generalizes to phase space.** Rather than "refuse when *time*
  curvature is large," the trustworthy-extrapolation test becomes "refuse when the
  `Ts(N)` relation changes character between the fit window and `N = 0`" — i.e. a
  detected kink in phase space, which *is* a feedback threshold (snowball
  bifurcation). This is a more physical gate than the time-domain curvature ratio.
- **This sharpens the ML note.** The learned target shifts from "infer τ" to
  "infer the phase-space feedback relation (its slope `λ`, and any curvature/kink)
  from a short burst + case descriptors" — a more physical and transferable target
  than a time constant.

## Open / next

- **A better nonlinear form.** The `1 − e^(N/N0)` form slightly overshoots past
  `N=0` on the concave case (pt10, Fig 12). A feedback-parameter model
  (`N = λ(T)·(T − T_eq)` with `λ` linear in `T`) or a two-parameter saturating
  form may tighten it and be more physically interpretable.
- **Wider validation.** Score against the other slow variables (energy_bot as a
  second imbalance coordinate; the ice reservoir) and with leave-one-run-out, per
  the ML note's discipline.
- **The residual-N endpoint bias.** Ground truth here is the final state at
  `N ≈ −1.3`, not true `N=0`; the ~0.15–0.3 K floor partly reflects that the runs
  never fully equilibrated. Real equilibria (or a longer reference run) would let
  us separate method error from endpoint bias.
- Everything above stays **offline**. The online run→jump→run harness (CLAUDE.md
  Open items O1/O2) remains out of scope; this note only establishes that the
  phase-space extrapolation is sound enough to eventually drive it.

---

## Ice growth is Stefan-limited (session 2026-09-25)

Trying to extrapolate the ice reservoir `hi` to `N = 0` (for the jump advisor)
failed badly: late in the runs the saturating fit wanted 700 m of ice. The
reason is physical, and it reframes what an ice "jump" is.

**The TOA deficit is being spent freezing ice, not leaking.** At year 148 the
surface imbalance `energy_bot` matches `energy_top` (both −0.6 to −2.1 W/m²),
and `hi` still grows 0.6–1.8 m/decade. Converting: 0.15 m/yr × ρ_i L
(≈ 3×10⁸ J/m³) ≈ 1.5 W/m² — the size of the residual N. The earlier docs'
"residual N ≈ −1.3 at the end of the runs" is not an endpoint bias; the runs
are still equilibrating, through the ice.

**`N × hi` is conserved.** In every case the product is constant to a few
percent from year 30 to 150 (e.g. pt13: −101.0, −100.8, −100.6, −100.3, −99.0
at yr 30/60/90/120/148). Fitting `N = a + b/hi` over yr 30–150 gives r² ≥ 0.98
in 12/15 cases (0.91–0.98 in the other three), with `a` ∈ [−0.05, +0.54] W/m².
This is conduction-limited (Stefan) growth: the deficit is conducted through
ice of thickness `h` (flux ∝ ΔT/h) and freezes onto its base, so `h²` grows
linearly in time and `N → a` only as `h → ∞`.

**Consequences.**

1. The ice never "equilibrates" at `N = 0` on any useful timescale — and
   `N = 0` is not a meaningful jump target for it. A jump instead *chooses a
   target imbalance* and sets the ice to the thickness the conduction law
   says produces it: `h_target = b / (N_target − a)`, ice factor
   `h_target / h_now`.
2. The right fit form for the reservoir is **hyperbolic**, not linear or
   saturating-exponential. Hindcast (truncate at year Y, 10-yr int2 window,
   predict `hi` at the future N the run actually reached): median capture of
   the true growth 0.90 to year 150, median error 2.0 m vs 6.2 m for linear
   and 10.8 m for not jumping.
3. It supplies a conversion from ice to *time*: since `h²` is linear in `t`,
   `years_skipped ≈ (h_new² − h_now²) / (dh²/dt)`.
4. TS stays well described by a linear `TS(N)` (Gregory); it is the
   *reference* the post-jump run is checked against, not a field that is
   written.
5. The Stefan law and the fixed-area `vicen` scaling both assume a settled ice
   edge. Windows ending at yr 20 (ICEFRAC still moving ≥ 0.03 per decade) were
   the only systematic failures (≈ 50 % error); by yr 30 |ΔICEFRAC| is ~0.002
   per decade. Hence the gate "|ΔICEFRAC| over the window ≤ 0.02".

pt10 (the one full snowball, ICEFRAC = 1, and the concave `Ts(N)` case in
Fig 12) has the weakest conduction-law fit and is refused by the advisor late
in its run.

## Advisor hindcast (defaults for `exocam-accelerate advise`)

Driver scripts (outside the repo): `../scratch/advisor_hindcast.py`,
`../scratch/advisor_endtoend.py`.

**Window.** Accuracy is nearly flat in window length (median relative `hi`
error 2–5 % for 10–40 yr windows at origins ≥ yr 30), but late in a run N moves
~0.1 W/m² per decade, so short windows lack N-range and fail `|corr| ≥ 0.9`
(yr 90–110 origins: 42 % accepted at 10 yr vs 84 % at 40 yr). Default: the
longest of 40/30/20/10 yr over which the ice edge has settled.

**End-to-end, default settings** (target = remove half the imbalance, ice
factor clipped at 1.5). Each run truncated at Y = 30…130; the advisor's
recommendation compared with the *real* run at the moment it actually reached
the recommended thickness:

| check | result |
|---|---|
| origins accepted | 144 / 165 (refusals concentrated at Y ≥ 100) |
| recommended factor | median 1.50 (the clip binds) |
| N_after − real N at arrival | median +0.09, p90 \|·\| 0.31 W/m² |
| TS predicted − real TS at arrival | median −0.09, p90 \|·\| 0.39 K |
| years_skipped / real years to arrival | median 0.95 (p10–p90 0.87–0.98) |

One capped jump ≈ 180 model years of ice growth for these cases.

**What this does not yet show.** It validates that the recommended state lies
on the trajectory the run actually follows. It does not show that *writing*
that ice into cice.r — with the SOM ocean, snow and atmosphere left as they
are — lets the model continue along it. That is the Tier-1/Tier-2 test
(restart-integration doc §6b), and the first real jumps are that test: after
the continuation, N should settle near the advice's `N_after` within a few
years.

## After a jump: advising again, and the post-jump check (2026-09-26)

**Advising again.** The int2 column is a 10-yr trailing mean, so for a decade
after a jump it straddles two states. `advise` detects a jump as a one-year
step in `hi` that stands out from its neighbours' growth (early spin-up grows
>8 %/yr for years and is not flagged), then fits only native annual means from
2 years after the jump. Hindcast (native annual means, hyperbolic fit, origins
yr 40–120): at |corr| ≥ 0.7 a 20-yr window gives median / p90 relative `hi`
error 2.7 % / 7.4 % — int2 accuracy — while 10-yr windows are noisy (7 % / 16 %).
Hence ≥ 15 settled post-jump years before the next jump. Late in a run the
window's N-range is small, so the target is clipped to 5x that range rather
than refused: later jumps come out smaller, the phase-space form of a
decreasing Δt schedule.

The ice-edge gate compares 3-yr means at the window's ends. A trend-based
version was tried and rejected: it passed windows reaching back into the
initial transient (yr ~10), which doubled the p90 `N_after` error
(1.09 vs 0.34 W/m²) in the end-to-end hindcast.

**The check** (`check.py`) scores the post-jump run against the advice's own
relations at the ice the run actually has, relative to the run's offset from
them over the 5 pre-jump years. Null-jump test (factor 1 at yr 40–130 on all
15 runs, `../scratch/check_falsealarm.py`): FAIL in 3/130 at 5 yr and 4/130 at
8 yr, all but one in pt10 — the full snowball, whose TS shifts by 1–2 K on its
own around yr 80–100 (flagging it is the conservative outcome). Sensitivity: a
jump the model rejects leaves N near its pre-jump value, N_now/3 away from
the law at a 1.5 factor, so a rejected jump is reliably caught when
|N_now| ≳ 1.2 W/m²; below that it is within the 0.4 W/m² noise floor.
