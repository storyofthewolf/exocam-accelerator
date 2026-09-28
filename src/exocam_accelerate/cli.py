"""Command line: ``exocam-accelerate {advise,jump,check,rollback,restore}``.

advise    read a case's exocam-trend .txt, fit the ice conduction law, print
          the recommended ice factor (--json saves it). Detects earlier jumps
          and then fits post-jump data only.
jump      pre-flight the run directory (no job queued/running, consistent
          rpointers, no output past the restart, pristine archived set), then
          scale vicen/eicen in the cice.r that rpointer.ice names, in place.
check     after the post-jump segment: PASS / WAIT / FAIL against the advice
          the jump was sized from (exit 0 / 10 / 20, for scripted polling).
rollback  reset the whole restart set to the pre-jump date (from the archive).
restore   before resubmitting only: put the pristine cice.r back.

Workflow per case: advise -> jump -> resubmit a SHORT continuation segment
(CONTINUE_RUN=TRUE) -> regenerate trends -> check -> PASS: resume normal
segments; FAIL: rollback and resubmit.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from .advise import AdvisorConfig, advise
from .phase_space import PhaseGateConfig
from .trend_io import case_start_year, load_case

EXIT_BLOCKED = 3


def _fmt(x, digits=4):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    return f"{x:.{digits}g}"


def _start_year(args) -> int:
    if args.start_year is not None:
        return args.start_year
    return case_start_year(args.trend_dir, args.case)


def cmd_advise(args) -> int:
    from .trend_io import file_provenance

    cols = load_case(args.trend_dir, args.case)
    provenance = file_provenance(args.trend_dir, args.case)
    gate = PhaseGateConfig(max_extrapolation_ratio=args.max_extrapolation_ratio)
    cfg = AdvisorConfig(window_years=args.window, which=args.which,
                        n_fraction=args.n_fraction, N_target=args.n_target,
                        max_ice_factor=args.max_ice_factor, gate=gate,
                        since_year=args.since)
    adv = advise(cols, args.case, cfg, start_year=_start_year(args),
                provenance=provenance)

    since = f"   post-jump since year {adv.since_year}" if adv.since_year else ""
    print(f"case {adv.case}   data through model year {adv.model_year}   fit "
          f"{adv.which_used}, last {adv.window_years:g} yr{since}")
    print(f"imbalance N: now {adv.N_now:+.2f} W/m2 (on conduction law "
          f"{_fmt(adv.N_now_fit, 3)})   target {_fmt(adv.N_target, 3)}")
    ice = adv.ice
    if ice is not None and "hyperbolic" in ice.fits:
        f = ice.fits["hyperbolic"]
        a, b = f.params
        print(f"ice law: N = {a:+.3f} {b:+.2f}/hi   corr(N,1/hi) {f.corr:+.3f}   "
              f"hi now {ice.X_now:.2f} m -> "
              f"{_fmt(ice.predictions.get('hyperbolic'))} m at target")
    print()
    for w in adv.warnings:
        print(f"warning: {w}")
    if not adv.jump:
        print("RECOMMENDATION: do not jump; keep running the model.")
        for r in adv.reasons:
            print(f"  - {r}")
    else:
        note = (f" (clipped from {adv.ice_factor_raw:.3f} at "
                f"{cfg.max_ice_factor:g})" if adv.ice_factor_clipped else "")
        print(f"RECOMMENDED ice factor (vicen, eicen): {adv.ice_factor:.4f}{note}")
        print(f"  hi {ice.X_now:.2f} -> {ice.X_now * adv.ice_factor:.2f} m; "
              f"skips ~{_fmt(adv.years_skipped, 3)} model years")
        print(f"  expected after adjustment: N ~ {adv.N_after:+.2f} W/m2", end="")
        for v, r in adv.temperatures.items():
            if r.accepted:
                print(f", {v} ~ {r.prediction:.2f} (now {r.X_now:.2f})", end="")
        print()

    if args.json:
        Path(args.json).write_text(json.dumps(adv.to_dict(), indent=2))
        print(f"wrote {args.json}")
    return 0


def _target(args) -> Path:
    from .restart import locate_cice_restart
    if args.cice_r:
        return Path(args.cice_r)
    return locate_cice_restart(args.rundir)


def _print_findings(findings) -> None:
    mark = {"ok": "  ok   ", "warn": "  WARN ", "block": "  BLOCK"}
    for f in findings:
        print(f"{mark[f.level]} {f.message}")


def cmd_jump(args) -> int:
    from .restart import apply_ice_jump
    from .runstate import active_jobs, blocked, preflight

    advice = None
    if args.advice:
        advice = json.loads(Path(args.advice).read_text())
        if advice.get("ice_factor") is None:
            print("advice says do not jump:", *advice.get("reasons", []), sep="\n  ")
            return 1
        ice_factor = float(advice["ice_factor"])
    else:
        ice_factor = args.ice_factor

    path = _target(args)
    probe = None if args.skip_slurm_check else active_jobs
    archive = Path(args.archive) if args.archive else None
    print("pre-flight:")
    findings = preflight(path, archive, probe, advice,
                         allow_no_archive=args.allow_no_archive_rollback)
    _print_findings(findings)
    if blocked(findings):
        print("refusing to jump: resolve the BLOCK items above")
        return EXIT_BLOCKED

    rec = apply_ice_jump(path, ice_factor, args.snow_factor, advice, dry_run=True)
    print(f"target   {rec.path}")
    print(f"backup   {rec.backup}"
          f"{' (exists; jump re-applied from it)' if rec.backup.exists() else ' (will be created)'}")
    print(f"factors  ice {rec.ice_factor:g}   snow {rec.snow_factor:g}   "
          f"first post-jump model year {rec.metadata['jump_model_year']}")
    for name, msg in rec.adjustments.items():
        print(f"  {name:>6}: {msg}")
    if args.dry_run:
        print("dry run: nothing written")
        return 0
    if not args.yes:
        if input("write this jump? [y/N] ").strip().lower() != "y":
            print("aborted")
            return 1
    rec = apply_ice_jump(path, ice_factor, args.snow_factor, advice)
    print(f"written and verified; jump log {rec.log}")
    if probe is not None:
        jobs = probe()
        if jobs and rec.metadata["case"] in jobs:
            print(f"WARNING: a job named {rec.metadata['case']} appeared while the "
                  f"jump was written — check it read the jumped file, or restore")
    print("next: resubmit a SHORT continuation segment (CONTINUE_RUN=TRUE), then "
          "regenerate trends and run 'exocam-accelerate check'")
    return 0


def cmd_check(args) -> int:
    from .check import CheckConfig, check_jump

    if args.log:
        log_file = Path(args.log)
    else:
        from .runstate import select_jump_log
        log_file = select_jump_log(Path(args.rundir).resolve(), args.date)
    log = json.loads(log_file.read_text())
    cols = load_case(args.trend_dir, args.case)
    cfg = CheckConfig(settle_years=args.settle_years, min_years=args.min_years)
    res = check_jump(cols, log, _start_year(args), cfg)

    print(f"case {args.case}   jump log {log_file.name}")
    print(f"jump at model year {res.jump_year} (factor {log['ice_factor']:g}); "
          f"{res.years_after} post-jump year(s), {res.settled_years} settled")
    m = res.metrics
    if "hi_first" in m:
        print(f"  landed: hi {m['hi_first']:.2f} m vs {m['hi_expected']:.2f} "
              f"expected ({m['land_error']:+.1%})")
    if "dN" in m:
        print(f"  N: observed {m['N_obs']:+.2f}, conduction law {m['N_law']:+.2f} "
              f"(diff {m['dN']:+.2f}, tolerance {m['tol_N']:.2f}) W/m2")
    if "dTS" in m:
        print(f"  TS: observed {m['TS_obs']:.2f} K (diff from Gregory relation "
              f"{m['dTS']:+.2f}, tolerance {m['tol_TS']:.2f})")
    if "hi_trend_m_per_yr" in m:
        print(f"  hi trend since jump {m['hi_trend_m_per_yr']:+.2f} m/yr")
    for r in res.reasons:
        print(f"  - {r}")
    advice_line = {"PASS": "jump holds: resume normal segments",
                   "WAIT": "keep the short segments going and check again",
                   "FAIL": "roll back: exocam-accelerate rollback --rundir ... --archive ..."}
    print(f"VERDICT: {res.verdict.value} — {advice_line[res.verdict.value]}")
    return res.verdict.exit_code


def cmd_rollback(args) -> int:
    from .runstate import active_jobs, blocked, execute_rollback, plan_rollback

    probe = None if args.skip_slurm_check else active_jobs
    archive = Path(args.archive) if args.archive else None
    plan = plan_rollback(args.rundir, archive, args.date, probe)
    print(f"rollback {plan.case} to {plan.date} (source: {plan.source})")
    _print_findings(plan.findings)
    for a in plan.actions:
        print(f"  will: {a}")
    if plan.discarded:
        print(f"  output from the discarded segment ({len(plan.discarded)} files) is "
              f"NOT deleted — move it aside before trend analysis, or it will be "
              f"overwritten as the rerun reaches those dates:")
        for d in plan.discarded[:10]:
            print(f"    {d}")
        if len(plan.discarded) > 10:
            print(f"    ... and {len(plan.discarded) - 10} more")
    if blocked(plan.findings):
        print("refusing to roll back: resolve the BLOCK items above")
        return EXIT_BLOCKED
    if args.dry_run:
        print("dry run: nothing changed")
        return 0
    if not args.yes:
        if input("roll back? [y/N] ").strip().lower() != "y":
            print("aborted")
            return 1
    execute_rollback(plan, archive)
    print(f"rolled back; resubmit as a continuation from {plan.date}")
    return 0


def cmd_restore(args) -> int:
    from .restart import restore
    path = _target(args)
    backup = restore(path)
    print(f"restored {path} from {backup}")
    return 0


def _add_target(p):
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--rundir", help="case run directory (cice.r found via rpointer.ice)")
    g.add_argument("--cice-r", help="explicit cice.r file")


def _add_trend(p):
    p.add_argument("trend_dir", help="directory holding <case>_*_{cam,cice}.txt")
    p.add_argument("case", help="case id (file-name prefix)")
    p.add_argument("--start-year", type=int, default=None,
                   help="model year of the first trend month (default: from the "
                        "trend file name)")


def _add_safety(p):
    p.add_argument("--archive", help="the case's short-term archive root "
                                     "(DOUT_S_ROOT, holding rest/<date>/); the "
                                     "rollback source, mandatory for 'jump' unless "
                                     "--allow-no-archive-rollback is also passed")
    p.add_argument("--skip-slurm-check", action="store_true",
                   help="skip the squeue probe (only after confirming yourself "
                        "that no job for the case is queued or running)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--yes", action="store_true", help="skip confirmation")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="exocam-accelerate", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("advise", help="recommend an ice factor from trend output")
    _add_trend(a)
    a.add_argument("--window", type=float, default=None,
                   help="fit window, years (default: longest of 40/30/20/10 "
                        "with a settled ice edge)")
    a.add_argument("--which", default="int2", choices=["native", "int1", "int2"],
                   help="trend column to fit before any jump (int2)")
    a.add_argument("--since", type=int, default=None,
                   help="first model year run from a jumped state (default: "
                        "detected from a step in hi)")
    a.add_argument("--n-fraction", type=float, default=0.5,
                   help="fraction of the current imbalance to remove (0.5)")
    a.add_argument("--n-target", type=float, default=None,
                   help="absolute target imbalance, W/m2 (overrides --n-fraction)")
    a.add_argument("--max-ice-factor", type=float, default=1.5,
                   help="hard clip on the ice factor (1.5)")
    a.add_argument("--max-extrapolation-ratio", type=float, default=5.0,
                   help="gate: max distance beyond the sampled N-range (5)")
    a.add_argument("--json", help="save advice as JSON (input to 'jump --advice')")
    a.set_defaults(func=cmd_advise)

    j = sub.add_parser("jump", help="pre-flight, then scale cice.r ice in place")
    _add_target(j)
    f = j.add_mutually_exclusive_group(required=True)
    f.add_argument("--advice", help="advice JSON from 'advise --json'")
    f.add_argument("--ice-factor", type=float,
                   help="explicit ice factor (no advice: 'check' cannot score it)")
    j.add_argument("--snow-factor", type=float, default=1.0,
                   help="scale vsnon/esnon too (default 1 = untouched)")
    _add_safety(j)
    j.add_argument("--allow-no-archive-rollback", action="store_true",
                   help="escape hatch for a missing --archive: proceed without an "
                        "archived restart set, after verifying that every "
                        "rpointer-named component restart is present in the run "
                        "directory (only a valid rollback source until the next "
                        "st_archive sweep — normally pass --archive instead)")
    j.set_defaults(func=cmd_jump)

    c = sub.add_parser("check", help="score the post-jump run: PASS/WAIT/FAIL")
    _add_trend(c)
    g = c.add_mutually_exclusive_group(required=True)
    g.add_argument("--rundir", help="run directory holding exocam_accelerate/*.accel.json")
    g.add_argument("--log", help="explicit jump log (.accel.json)")
    c.add_argument("--date", help="restart date, if several jumps are active")
    c.add_argument("--settle-years", type=int, default=2)
    c.add_argument("--min-years", type=int, default=3)
    c.set_defaults(func=cmd_check)

    r = sub.add_parser("rollback", help="reset the restart set to the pre-jump date")
    r.add_argument("--rundir", required=True)
    r.add_argument("--date", help="restart date, if several jumps are active")
    _add_safety(r)
    r.set_defaults(func=cmd_rollback)

    s = sub.add_parser("restore", help="before resubmitting: pristine cice.r back")
    _add_target(s)
    s.set_defaults(func=cmd_restore)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (FileNotFoundError, KeyError, ValueError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
