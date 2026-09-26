"""Command line: ``exocam-accelerate {advise,jump,restore}``.

advise   read a case's exocam-trend .txt, print the phase-space extrapolation
         per slow variable and the recommended ice factor (--json saves it).
jump     apply an ice factor (from --advice JSON or --ice-factor) to the
         cice.r that rpointer.ice names, in place with a pristine backup.
restore  put the pristine cice.r back.

After a jump, resubmit the case as a continuation (CONTINUE_RUN=TRUE).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from .advise import AdvisorConfig, advise
from .phase_space import PhaseGateConfig
from .trend_io import load_case


def _fmt(x, digits=4):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    return f"{x:.{digits}g}"


def cmd_advise(args) -> int:
    cols = load_case(args.trend_dir, args.case)
    gate = PhaseGateConfig(max_extrapolation_ratio=args.max_extrapolation_ratio)
    cfg = AdvisorConfig(window_years=args.window, which=args.which,
                        n_fraction=args.n_fraction, N_target=args.n_target,
                        max_ice_factor=args.max_ice_factor, gate=gate)
    adv = advise(cols, args.case, cfg)

    print(f"case {adv.case}   year {adv.year_now:.0f}   fit {cfg.which}, "
          f"last {adv.window_years:g} yr")
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


def cmd_jump(args) -> int:
    from .restart import apply_ice_jump

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
    if advice and advice.get("case") and advice["case"] not in path.name:
        print(f"warning: advice is for case {advice['case']!r} but target is "
              f"{path.name}")

    rec = apply_ice_jump(path, ice_factor, args.snow_factor, advice, dry_run=True)
    print(f"target   {rec.path}")
    print(f"backup   {rec.backup.name}"
          f"{' (exists; jump re-applied from it)' if rec.backup.exists() else ' (will be created)'}")
    print(f"factors  ice {rec.ice_factor:g}   snow {rec.snow_factor:g}")
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
    print(f"written and verified; log {rec.path.name}.accel.json")
    print("next: resubmit the case as a continuation (CONTINUE_RUN=TRUE)")
    return 0


def cmd_restore(args) -> int:
    from .restart import restore
    path = _target(args)
    backup = restore(path)
    print(f"restored {path} from {backup.name}")
    return 0


def _add_target(p):
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--rundir", help="case run directory (cice.r found via rpointer.ice)")
    g.add_argument("--cice-r", help="explicit cice.r file")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="exocam-accelerate", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("advise", help="recommend an ice factor from trend output")
    a.add_argument("trend_dir", help="directory holding <case>_*_{cam,cice}.txt")
    a.add_argument("case", help="case id (file-name prefix)")
    a.add_argument("--window", type=float, default=None,
                   help="fit window, years (default: longest of 40/30/20/10 "
                        "with a settled ice edge)")
    a.add_argument("--which", default="int2", choices=["native", "int1", "int2"],
                   help="trend column to fit (int2)")
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

    j = sub.add_parser("jump", help="apply an ice factor to cice.r in place")
    _add_target(j)
    f = j.add_mutually_exclusive_group(required=True)
    f.add_argument("--advice", help="advice JSON from 'advise --json'")
    f.add_argument("--ice-factor", type=float, help="explicit ice factor")
    j.add_argument("--snow-factor", type=float, default=1.0,
                   help="scale vsnon/esnon too (default 1 = untouched)")
    j.add_argument("--dry-run", action="store_true")
    j.add_argument("--yes", action="store_true", help="skip confirmation")
    j.set_defaults(func=cmd_jump)

    r = sub.add_parser("restore", help="restore the pristine cice.r")
    _add_target(r)
    r.set_defaults(func=cmd_restore)
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
