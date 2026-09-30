"""Command line: ``exocam-accelerate {advise,taper,advise-ocean,pattern,jump,check,
rollback,restore,somtp-map,view}``.

Cold (sea-ice) regime — aqua_ice plugin, cice.r:

advise    read a case's exocam-trend .txt, fit the ice conduction law, print
          the recommended ice factor (--json saves it). Detects earlier jumps
          and then fits post-jump data only.
taper     re-size that advice for a tapered jump: from two archived restarts,
          weight each cell by how conduction-limited its ice growth is, and
          solve for the peak factor that still reaches the target imbalance
          (writes tapered advice JSON + a .taper.nc weight map; needs netCDF4).
jump      pre-flight the run directory (no job queued/running, consistent
          rpointers, no output past the restart, pristine archived set), then
          scale vicen/eicen in the cice.r that rpointer.ice names, in place.
check     after the post-jump segment: PASS / WAIT / FAIL against the advice
          the jump was sized from (exit 0 / 10 / 20, for scripted polling).
rollback  reset the whole restart set to the pre-jump date (from the archive).
restore   before resubmitting only: put the pristine cice.r back.
view      local interactive viewer (http://127.0.0.1:8765) of the trend files and
          jump logs in a directory: phase space, fits, jumps, check verdicts.

Hot, ice-free regime — som_ocean plugin, docn.r somtp:
advise-ocean  fit TS against the surface imbalance energy_bot (the slab ocean's
              own equilibrium; energy_top also carries the heat the atmosphere
              stores while it warms), print the recommended somtp increment,
              scaled by the measured heat ratio C_total/C_ocean
              print the recommended somtp increment (--json saves it).
pattern       optional: shape the increment by the local warming rate between
              two archived docn.r restarts (area mean unchanged).
somtp-map     write docn.r somtp as a lat-lon netCDF (ncview / Panoply / view).
jump / check / rollback / restore work on either plugin; jump picks the file
from the advice (rpointer.ocn for ocean advice).

Workflow per case: advise[-ocean] -> [taper | pattern] -> jump -> resubmit a
SHORT continuation segment (CONTINUE_RUN=TRUE) -> regenerate trends -> check
-> PASS: resume normal segments; FAIL: rollback and resubmit.
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


def cmd_taper(args) -> int:
    from .advise import taper_advice
    from .restart import build_taper_mask, first_model_year, write_taper_file
    from .taper import TaperConfig

    advice = json.loads(Path(args.advice).read_text())
    if advice.get("ice_factor") is None:
        print("advice says do not jump:", *advice.get("reasons", []), sep="\n  ")
        return 1
    case = advice["case"]
    date = args.date or f"{int(advice['model_year']) + 1:04d}-01-01-00000"
    if first_model_year(date) - 1 != int(advice["model_year"]):
        print(f"warning: advice data run through model year {advice['model_year']}, "
              f"the taper restart is {date}")
    cfg = TaperConfig(ramp_lo=args.ramp[0], ramp_hi=args.ramp[1])
    mask, sources = build_taper_mask(args.archive, case, date, args.baseline_years,
                                     args.grid_file, cfg)
    out = taper_advice(advice, mask, args.max_ice_factor)

    out_json = Path(args.json)
    stem = out_json.name[:-5] if out_json.name.endswith(".json") else out_json.name
    weights = out_json.with_name(stem + ".taper.nc")
    sha = write_taper_file(weights, mask, out["ice_factor"])
    out["taper"].update(restart_date=date, weights_file=str(weights.resolve()),
                        weights_sha256=sha, sources=sources)

    t, u = out["taper"], out["taper"]["uniform"]
    print(f"case {case}   restart {date}   baseline {args.baseline_years} yr   "
          f"S_ref {t['S_ref_m2_per_yr']:.2f} m2/yr   ramp {cfg.ramp_lo:g}-{cfg.ramp_hi:g}")
    print(f"area: full weight {t['area_full_weight']:.1%}   partial "
          f"{t['area_partial_weight']:.1%}   unscaled {t['area_unscaled']:.1%}")
    note = (f" (clipped from {_fmt(out['ice_factor_raw'], 4)})"
            if out["ice_factor_clipped"] else "")
    print(f"peak factor {out['ice_factor']:.4f}{note}   global-mean hi x"
          f"{out['effective_factor']:.4f}")
    print(f"expected N after {out['N_after']:+.2f} W/m2 (uniform x{u['ice_factor']:g}: "
          f"{u['N_after']:+.2f})   skips ~{_fmt(out['years_skipped'], 3)} yr "
          f"(uniform ~{_fmt(u['years_skipped'], 3)})")
    for v, x in (out.get("expected_after_jump") or {}).items():
        if x is not None:
            print(f"  {v} ~ {x:.2f}")
    for w in out["warnings"][len(advice.get("warnings") or []):]:
        print(f"warning: {w}")
    out_json.write_text(json.dumps(out, indent=2))
    print(f"wrote {out_json} and {weights}")
    return 0


def _since_from_logs(rundir, plugin: str):
    """Latest first post-jump model year among a run directory's active jump
    logs for ``plugin`` (None when there are none)."""
    from .restart import find_jump_logs
    years = []
    for p in find_jump_logs(rundir):
        log = json.loads(p.read_text())
        if log.get("plugin") == plugin:
            years.append(int(log["jump_model_year"]))
    return max(years) if years else None


def _print_probe(d, json_path) -> int:
    print(f"case {d['case']}   data through model year {d['model_year']}   PROBE mode "
          f"(maps the trajectory; not sized to an equilibrium)")
    if "TS_now" in d:
        print(f"TS now {d['TS_now']:.2f} K, trend {d['trend_K_per_yr']:+.3f} K/yr "
              f"(sigma {d['sigma_TS']:.2f} K)   {d['config']['imbalance']} "
              f"{d['N_now']:+.2f} W/m2 (last {d['config']['recent_years']} yr, "
              f"sigma {d['sigma_N']:.2f})")
    lr = d.get("lambda_recent")
    if lr:
        print(f"recent lambda (two {d['config']['recent_years']}-yr periods, "
              f"dTS {lr['dTS']:+.2f} K): {lr['value']:+.2f} ± {lr['se']:.2f} W/m2/K")
    print()
    for w in d["warnings"]:
        print(f"warning: {w}")
    if d["somtp_dT"] is None:
        print("RECOMMENDATION: no probe.")
        for r in d["reasons"]:
            print(f"  - {r}")
    else:
        how = (f"{d['config']['probe_years']:g} yr of the recent trend"
               if d.get("sizing") == "trend" else "explicit --probe-dt")
        print(f"PROBE somtp increment: {d['somtp_dT']:+.3f} K ({how}; heat ratio "
              f"{d['heat']['used']:.2f})")
        print(f"  TS {d['TS_now']:.2f} -> {d['TS_after']:.2f} K ({d['dTS']:+.2f}); after "
              f"{d['config']['recent_years']} settled years 'check' resolves "
              f"lambda >= ~{d['lambda_detectable']:.2f} W/m2/K")
    if json_path:
        Path(json_path).write_text(json.dumps(d, indent=2))
        print(f"wrote {json_path}")
    return 0


def cmd_advise_ocean(args) -> int:
    from .ocean_advise import OceanAdvisorConfig, ProbeConfig, advise_ocean, probe_ocean
    from .trend_io import file_provenance

    cols = load_case(args.trend_dir, args.case)
    since = args.since
    if since is None and args.rundir:
        since = _since_from_logs(args.rundir, "som_ocean")
        if since is not None:
            print(f"jump log in {args.rundir}: fitting post-jump data since model "
                  f"year {since}")
    if args.probe or args.probe_dt is not None:
        pc = ProbeConfig(probe_years=args.probe_years, probe_dT=args.probe_dt,
                         max_dT=args.max_dt, imbalance=args.imbalance, since_year=since,
                         heat_ratio=args.heat_ratio)
        d = probe_ocean(cols, args.case, pc, start_year=_start_year(args),
                        provenance=file_provenance(args.trend_dir, args.case))
        return _print_probe(d, args.json)
    gate = PhaseGateConfig(max_extrapolation_ratio=args.max_extrapolation_ratio)
    cfg = OceanAdvisorConfig(window_years=args.window, which=args.which,
                             imbalance=args.imbalance, n_fraction=args.n_fraction,
                             N_target=args.n_target, max_dT=args.max_dt,
                             heat_ratio=args.heat_ratio, gate=gate, since_year=since)
    adv = advise_ocean(cols, args.case, cfg, start_year=_start_year(args),
                       provenance=file_provenance(args.trend_dir, args.case))
    d = adv.to_dict()
    since_s = f"   post-jump since year {adv.since_year}" if adv.since_year else ""
    print(f"case {adv.case}   data through model year {adv.model_year}   fit "
          f"{adv.which_used}, last {adv.window_years:g} yr{since_s}")
    print(f"TS now {adv.TS_now:.2f} K   {cfg.imbalance}: now {adv.N_now:+.2f} W/m2 "
          f"(on the line {_fmt(adv.N_now_fit, 3)})   target {_fmt(adv.N_target, 3)}")
    g = adv.gregory
    if g is not None and "linear" in g.fits:
        c0, c1 = g.fits["linear"].params
        print(f"Gregory line: TS = {c0:.2f} {c1:+.3f}*N   corr {g.fits['linear'].corr:+.3f}"
              f"   TS at N=0: {c0:.2f} K")
    lk = d["gap"]
    if lk.get("mean") is not None:
        print(f"energy_top - energy_bot over the window: {lk['mean']:+.2f} W/m2 "
              f"(heat stored by the atmosphere; trend {lk['trend_per_yr']:+.3f}/yr)")
    h = d["heat"]
    if h.get("measured"):
        print(f"heat capacities: ocean {h['C_ocean']:.1f}, total {h['C_total']:.1f} "
              f"W yr/m2/K -> heat ratio {h['ratio']:.2f} (used {h['used']:.2f})")
    if d["C_eff_m_seawater"] is not None:
        print(f"effective heat capacity {d['C_eff_W_yr_m2_K']:.2f} W yr/m2/K "
              f"(~{d['C_eff_m_seawater']:.0f} m of sea water)   relaxation time "
              f"~{_fmt(adv.tau_years, 3)} yr")
    print()
    for w in adv.warnings:
        print(f"warning: {w}")
    if not adv.jump:
        print("RECOMMENDATION: do not jump; keep running the model.")
        for r in adv.reasons:
            print(f"  - {r}")
    else:
        note = (f" (clipped from {adv.somtp_dT_raw:+.3f})" if adv.clipped else "")
        print(f"RECOMMENDED somtp increment: {adv.somtp_dT:+.3f} K{note}"
              f"   (heat ratio {adv.heat['used']:.2f})")
        print(f"  TS {adv.TS_now:.2f} -> {adv.TS_after:.2f} K; skips ~"
              f"{_fmt(adv.years_skipped, 3)} model years")
        print(f"  expected after adjustment: {cfg.imbalance} ~ {adv.N_after:+.2f} W/m2")
    if args.json:
        Path(args.json).write_text(json.dumps(d, indent=2))
        print(f"wrote {args.json}")
    return 0


def cmd_pattern(args) -> int:
    from .ocean_advise import pattern_advice
    from .restart import build_ocean_pattern, first_model_year, write_pattern_file
    from .som_ocean import PatternConfig

    advice = json.loads(Path(args.advice).read_text())
    if advice.get("somtp_dT") is None:
        print("advice says do not jump:", *advice.get("reasons", []), sep="\n  ")
        return 1
    case = advice["case"]
    date = args.date or f"{int(advice['model_year']) + 1:04d}-01-01-00000"
    if first_model_year(date) - 1 != int(advice["model_year"]):
        print(f"warning: advice data run through model year {advice['model_year']}, "
              f"the pattern restart is {date}")
    cfg = PatternConfig(smooth=args.smooth, max_weight=args.max_weight)
    pat, grid, sources = build_ocean_pattern(args.archive, case, date, args.domain_file,
                                             args.baseline_years, cfg)
    summary = pat.summary()
    out_json = Path(args.json)
    stem = out_json.name[:-5] if out_json.name.endswith(".json") else out_json.name
    weights = out_json.with_name(stem + ".pattern.nc")
    sha = write_pattern_file(weights, pat, grid, advice["somtp_dT"])
    summary.update(restart_date=date, weights_file=str(weights.resolve()),
                   weights_sha256=sha, sources=sources,
                   domain_file=str(Path(args.domain_file).resolve()))
    out = pattern_advice(advice, summary)
    print(f"case {case}   restart {date}   baseline {args.baseline_years} yr   "
          f"mean somtp rate {summary['mean_rate_K_per_yr']:+.3f} K/yr")
    print(f"weights {summary['weight_min']:.2f}-{summary['weight_max']:.2f} "
          f"(std {summary['weight_std']:.2f}; area at the bounds "
          f"{summary['area_at_min']:.1%} / {summary['area_at_max']:.1%})")
    print(f"per-cell increment {advice['somtp_dT'] * summary['weight_min']:+.2f} to "
          f"{advice['somtp_dT'] * summary['weight_max']:+.2f} K, area mean "
          f"{advice['somtp_dT']:+.2f} K")
    out_json.write_text(json.dumps(out, indent=2))
    print(f"wrote {out_json} and {weights}")
    return 0


def cmd_somtp_map(args) -> int:
    from .restart import find_docn_domain, locate_docn_restart, write_somtp_map
    src = Path(args.docn_r) if args.docn_r else locate_docn_restart(args.rundir)
    domain = args.domain_file or find_docn_domain(src.parent)
    if domain is None:
        raise ValueError("no --domain-file, and no readable domainfile in the run "
                         "directory's docn_ocn_in")
    out = Path(args.out) if args.out else Path(src.name[:-3] + ".latlon.nc")
    write_somtp_map(src, domain, out)
    print(f"wrote {out} (somtp(lat, lon) in K from {src.name})")
    return 0


def _is_ocean(advice) -> bool:
    return bool(advice) and advice.get("plugin") == "som_ocean"


def _target(args, ocean: bool = False) -> Path:
    from .restart import locate_cice_restart, locate_docn_restart
    if getattr(args, "cice_r", None):
        return Path(args.cice_r)
    if getattr(args, "docn_r", None):
        return Path(args.docn_r)
    return locate_docn_restart(args.rundir) if ocean else locate_cice_restart(args.rundir)


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
    ocean = _is_ocean(advice) or args.delta_t is not None or bool(args.docn_r)
    if ocean:
        return _jump_ocean(args, advice)
    if advice is not None:
        if advice.get("ice_factor") is None:
            print("advice says do not jump:", *advice.get("reasons", []), sep="\n  ")
            return 1
        ice_factor = float(advice["ice_factor"])
    else:
        ice_factor = args.ice_factor
        if ice_factor is None:
            raise ValueError("give --advice, --ice-factor or --delta-t")

    path = _target(args)
    taper_file = None
    if advice is not None and advice.get("taper"):
        from .restart import file_sha256, restart_date
        t = advice["taper"]
        taper_file = Path(args.taper_file or t.get("weights_file") or "")
        if not taper_file.is_file():
            raise FileNotFoundError(f"tapered advice: weight map {taper_file} not found "
                                    f"(pass --taper-file)")
        if file_sha256(taper_file) != t.get("weights_sha256"):
            raise RuntimeError(f"{taper_file.name} does not match the sha256 recorded "
                               f"in the advice — refusing")
        if t.get("restart_date") != restart_date(path):
            raise RuntimeError(f"tapered advice was built for restart "
                               f"{t.get('restart_date')}, target is {restart_date(path)}")
    probe = None if args.skip_slurm_check else active_jobs
    archive = Path(args.archive) if args.archive else None
    print("pre-flight:")
    findings = preflight(path, archive, probe, advice,
                         allow_no_archive=args.allow_no_archive_rollback)
    _print_findings(findings)
    if blocked(findings):
        print("refusing to jump: resolve the BLOCK items above")
        return EXIT_BLOCKED

    rec = apply_ice_jump(path, ice_factor, args.snow_factor, advice, dry_run=True,
                         taper_file=taper_file)
    print(f"target   {rec.path}")
    print(f"backup   {rec.backup}"
          f"{' (exists; jump re-applied from it)' if rec.backup.exists() else ' (will be created)'}")
    print(f"factors  ice {rec.ice_factor:g}   snow {rec.snow_factor:g}   "
          f"first post-jump model year {rec.metadata['jump_model_year']}")
    if taper_file is not None:
        tm = rec.metadata["taper"]
        print(f"tapered  cell factors {tm['cell_factor_min']:.3f}-"
              f"{tm['cell_factor_max']:.3f}   global-mean hi "
              f"x{advice['effective_factor']:.4f}   map {taper_file}")
    for name, msg in rec.adjustments.items():
        print(f"  {name:>6}: {msg}")
    if args.dry_run:
        print("dry run: nothing written")
        return 0
    if not args.yes:
        if input("write this jump? [y/N] ").strip().lower() != "y":
            print("aborted")
            return 1
    rec = apply_ice_jump(path, ice_factor, args.snow_factor, advice,
                         taper_file=taper_file)
    print(f"written and verified; jump log {rec.log}")
    if probe is not None:
        jobs = probe()
        if jobs and rec.metadata["case"] in jobs:
            print(f"WARNING: a job named {rec.metadata['case']} appeared while the "
                  f"jump was written — check it read the jumped file, or restore")
    print("next: resubmit a SHORT continuation segment (CONTINUE_RUN=TRUE), then "
          "regenerate trends and run 'exocam-accelerate check'")
    return 0


def _jump_ocean(args, advice) -> int:
    from .restart import apply_ocean_jump, file_sha256, restart_date
    from .runstate import active_jobs, blocked, preflight

    if args.cice_r or args.ice_factor is not None or args.taper_file:
        raise ValueError("ocean jump: --cice-r / --ice-factor / --taper-file do not apply")
    if advice is not None:
        if not _is_ocean(advice):
            raise ValueError("--docn-r / --delta-t given with ice advice")
        if advice.get("somtp_dT") is None:
            print("advice says do not jump:", *advice.get("reasons", []), sep="\n  ")
            return 1
        dT = float(advice["somtp_dT"])
    else:
        dT = args.delta_t
    from .restart import find_docn_domain
    path = _target(args, ocean=True)
    pattern_file = None
    domain_file = args.domain_file or find_docn_domain(path.parent)
    if advice is not None and advice.get("pattern"):
        pt = advice["pattern"]
        pattern_file = Path(args.pattern_file or pt.get("weights_file") or "")
        if not pattern_file.is_file():
            raise FileNotFoundError(f"patterned advice: weight map {pattern_file} not "
                                    f"found (pass --pattern-file)")
        if file_sha256(pattern_file) != pt.get("weights_sha256"):
            raise RuntimeError(f"{pattern_file.name} does not match the sha256 recorded "
                               f"in the advice — refusing")
        if pt.get("restart_date") != restart_date(path):
            raise RuntimeError(f"patterned advice was built for restart "
                               f"{pt.get('restart_date')}, target is {restart_date(path)}")
        domain_file = domain_file or pt.get("domain_file")
    elif args.pattern_file:
        raise ValueError("--pattern-file needs patterned advice (exocam-accelerate pattern)")

    probe = None if args.skip_slurm_check else active_jobs
    archive = Path(args.archive) if args.archive else None
    print("pre-flight:")
    findings = preflight(path, archive, probe, advice,
                         allow_no_archive=args.allow_no_archive_rollback)
    _print_findings(findings)
    if blocked(findings):
        print("refusing to jump: resolve the BLOCK items above")
        return EXIT_BLOCKED
    rec = apply_ocean_jump(path, dT, advice, dry_run=True, pattern_file=pattern_file,
                           domain_file=domain_file)
    m = rec.metadata
    print(f"target   {rec.path}")
    print(f"backup   {rec.backup}"
          f"{' (exists; jump re-applied from it)' if rec.backup.exists() else ' (will be created)'}")
    print(f"somtp    {dT:+.3f} K requested, open-ocean mean {m['somtp_dT_applied_mean']:+.3f} K "
          f"(cells {m['cell_dT_min']:+.2f}..{m['cell_dT_max']:+.2f})   first post-jump "
          f"model year {m['jump_model_year']}")
    if domain_file is None:
        print("         no --domain-file (and none in docn_ocn_in): unweighted mean, "
              "every cell treated as ocean")
    else:
        print(f"domain   {domain_file}")
    for name, msg in rec.adjustments.items():
        print(f"  {name:>6}: {msg}")
    if args.dry_run:
        print("dry run: nothing written")
        return 0
    if not args.yes:
        if input("write this jump? [y/N] ").strip().lower() != "y":
            print("aborted")
            return 1
    rec = apply_ocean_jump(path, dT, advice, pattern_file=pattern_file,
                           domain_file=domain_file)
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
    from .check import check_any

    if args.log:
        log_file = Path(args.log)
    else:
        from .runstate import select_jump_log
        log_file = select_jump_log(Path(args.rundir).resolve(), args.date)
    log = json.loads(log_file.read_text())
    cols = load_case(args.trend_dir, args.case)
    res = check_any(cols, log, _start_year(args), args.settle_years, args.min_years)

    print(f"case {args.case}   jump log {log_file.name}")
    what = (f"somtp {log['somtp_dT']:+g} K" if log.get("plugin") == "som_ocean"
            else f"factor {log['ice_factor']:g}")
    print(f"jump at model year {res.jump_year} ({what}); "
          f"{res.years_after} post-jump year(s), {res.settled_years} settled")
    m = res.metrics
    if "dTS_first" in m:
        print(f"  landed: first-year TS moved {m['dTS_first']:+.2f} K net of drift vs "
              f"{m['dTS_expected']:+.2f} expected"
              + (f" (implied heat ratio {m['heat_ratio_implied']:.2f})"
                 if "heat_ratio_implied" in m else ""))
    if "lambda" in m:
        print(f"  probe: TS {m['TS_pre']:.2f} -> {m['TS_post']:.2f} K, energy "
              f"{m['N_pre']:+.2f} -> {m['N_post']:+.2f} W/m2 (± {m['se_dN']:.2f}); "
              f"lambda {m['lambda']:+.2f} ± {m['lambda_se']:.2f} W/m2/K"
              + (f"; TS_eq ~{m['TS_eq']:.1f} K" if "TS_eq" in m else ""))
    if "N_line" in m:
        print(f"  N: observed {m['N_obs']:+.2f}, Gregory line {m['N_line']:+.2f} "
              f"(diff {m['dN']:+.2f}, tolerance {m['tol_N']:.2f}) W/m2; TS "
              f"{m['TS_obs']:.2f} K")
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
    path = _target(args, ocean=args.ocean)
    backup = restore(path)
    print(f"restored {path} from {backup}")
    return 0


def _add_target(p):
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--rundir", help="case run directory (cice.r found via rpointer.ice, "
                                    "docn.r via rpointer.ocn for ocean advice)")
    g.add_argument("--cice-r", help="explicit cice.r file")
    g.add_argument("--docn-r", help="explicit docn.r file (ocean jump)")


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


def cmd_view(args) -> int:
    from .viewer import serve

    serve(args.directory, host=args.host, port=args.port,
          open_browser=not args.no_browser)
    return 0


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

    t = sub.add_parser("taper", help="re-size advice for a tapered (per-cell) jump")
    t.add_argument("--advice", required=True, help="advice JSON from 'advise --json'")
    t.add_argument("--archive", required=True,
                   help="the case's short-term archive root (rest/<date>/, ice/hist/)")
    t.add_argument("--date", help="restart date to jump (default: the January after "
                                  "the advice's last model year)")
    t.add_argument("--baseline-years", type=int, default=10,
                   help="years between the two archived restarts the growth "
                        "rates come from (10)")
    t.add_argument("--ramp", type=float, nargs=2, default=(0.05, 0.40),
                   metavar=("LO", "HI"),
                   help="weight 0 below LO, 1 above HI, in S/S_ref (0.05 0.40)")
    t.add_argument("--max-ice-factor", type=float, default=None,
                   help="clip on the peak factor (default: the advice's)")
    t.add_argument("--grid-file", help="cice.h file for tarea/tmask (default: "
                                       "latest in <archive>/ice/hist)")
    t.add_argument("--json", required=True,
                   help="tapered advice JSON to write (the weight map goes "
                        "beside it as <name>.taper.nc)")
    t.set_defaults(func=cmd_taper)

    o = sub.add_parser("advise-ocean",
                       help="recommend a docn.r somtp increment (hot, ice-free runs)")
    _add_trend(o)
    o.add_argument("--window", type=float, default=None,
                   help="fit window, years (default: longest of 40/30/20/10 whose "
                        "TS(N) line passes the gate and the current state)")
    o.add_argument("--which", default="int2", choices=["native", "int1", "int2"])
    o.add_argument("--imbalance", default="energy_bot",
                   choices=["energy_bot", "energy_top"],
                   help="phase-space coordinate (energy_bot: the slab's own "
                        "equilibrium; energy_top includes atmospheric storage)")
    o.add_argument("--since", type=int, default=None,
                   help="first model year run from a jumped state")
    o.add_argument("--rundir", help="read --since from the active ocean jump log here")
    o.add_argument("--n-fraction", type=float, default=0.5,
                   help="fraction of the current imbalance to remove (0.5)")
    o.add_argument("--n-target", type=float, default=None,
                   help="absolute target imbalance, W/m2 (overrides --n-fraction)")
    o.add_argument("--max-dt", type=float, default=10.0,
                   help="hard clip on the somtp increment, K (10)")
    o.add_argument("--heat-ratio", type=float, default=None,
                   help="somtp increment / TS change, to cover the heat the "
                        "atmosphere takes back (default: measured C_total/C_ocean, "
                        "clipped to [1, 4]; 'check' reports the implied value)")
    o.add_argument("--max-extrapolation-ratio", type=float, default=5.0)
    o.add_argument("--probe", action="store_true",
                   help="probe mode: step TS ahead by the recent trend x "
                        "--probe-years to map the trajectory (when the Gregory "
                        "line is not constrained)")
    o.add_argument("--probe-years", type=float, default=15.0)
    o.add_argument("--probe-dt", type=float, default=None,
                   help="explicit probe size, K (implies --probe; skips the trend "
                        "gate)")
    o.add_argument("--json", help="save advice as JSON (input to 'jump --advice')")
    o.set_defaults(func=cmd_advise_ocean)

    pt = sub.add_parser("pattern", help="shape an ocean jump by the measured warming "
                                        "pattern (optional)")
    pt.add_argument("--advice", required=True, help="advice JSON from 'advise-ocean'")
    pt.add_argument("--archive", required=True, help="short-term archive root (rest/)")
    pt.add_argument("--domain-file", required=True,
                    help="docn domain (or pop_frc) file: mask, area, xc, yc")
    pt.add_argument("--date", help="restart date (default: the January after the "
                                   "advice's last model year)")
    pt.add_argument("--baseline-years", type=int, default=10)
    pt.add_argument("--smooth", type=int, default=1,
                    help="box-smoothing half-width, cells (1 = 3x3)")
    pt.add_argument("--max-weight", type=float, default=3.0)
    pt.add_argument("--json", required=True,
                    help="patterned advice JSON (map beside it as <name>.pattern.nc)")
    pt.set_defaults(func=cmd_pattern)

    m = sub.add_parser("somtp-map", help="docn.r somtp -> lat-lon netCDF")
    mg = m.add_mutually_exclusive_group(required=True)
    mg.add_argument("--rundir", help="run directory (docn.r via rpointer.ocn)")
    mg.add_argument("--docn-r", help="explicit docn.r file")
    m.add_argument("--domain-file",
                   help="docn domain (or pop_frc) file: mask, area, xc, yc "
                        "(default: domainfile in the run directory's docn_ocn_in)")
    m.add_argument("-o", "--out", help="output (default <docn.r name>.latlon.nc here)")
    m.set_defaults(func=cmd_somtp_map)

    j = sub.add_parser("jump", help="pre-flight, then edit cice.r ice or docn.r "
                                    "somtp in place")
    _add_target(j)
    f = j.add_mutually_exclusive_group(required=True)
    f.add_argument("--advice", help="advice JSON from 'advise' / 'advise-ocean' "
                                    "(or taper / pattern)")
    f.add_argument("--ice-factor", type=float,
                   help="explicit ice factor (no advice: 'check' cannot score it)")
    f.add_argument("--delta-t", type=float,
                   help="explicit somtp increment, K (no advice: 'check' cannot "
                        "score it)")
    j.add_argument("--domain-file", help="ocean jump: docn domain file (mask, area; "
                                         "default: domainfile in docn_ocn_in)")
    j.add_argument("--pattern-file", help="weight map for patterned ocean advice "
                                          "(default: the path recorded in it)")
    j.add_argument("--taper-file", help="weight map for tapered advice (default: "
                                        "the path recorded in the advice)")
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
    c.add_argument("--settle-years", type=int, default=None, help="(2)")
    c.add_argument("--min-years", type=int, default=None, help="(3)")
    c.set_defaults(func=cmd_check)

    r = sub.add_parser("rollback", help="reset the restart set to the pre-jump date")
    r.add_argument("--rundir", required=True)
    r.add_argument("--date", help="restart date, if several jumps are active")
    _add_safety(r)
    r.set_defaults(func=cmd_rollback)

    s = sub.add_parser("restore", help="before resubmitting: pristine cice.r / "
                                       "docn.r back")
    _add_target(s)
    s.add_argument("--ocean", action="store_true",
                   help="with --rundir: restore the docn.r (rpointer.ocn)")
    s.set_defaults(func=cmd_restore)

    v = sub.add_parser("view", help="interactive local viewer of trends, fits and jumps")
    v.add_argument("directory", help="directory holding <case>_*_{cam,cice}.txt and "
                                     "any *.accel.json jump logs (subdirs searched too)")
    v.add_argument("--port", type=int, default=8765)
    v.add_argument("--host", default="127.0.0.1",
                   help="bind address (default localhost only)")
    v.add_argument("--no-browser", action="store_true", help="do not open a browser")
    v.set_defaults(func=cmd_view)
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
