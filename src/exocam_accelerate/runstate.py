"""Run-directory safety: pre-flight checks before a jump, and whole-set rollback.

No netCDF needed except to confirm an archived cice.r is pristine (optional).

Facts about CESM 1.2.1 this relies on (scripts/ccsm_utils/Tools/st_archive.sh):
* At the end of every successful segment st_archive moves the newest restart
  of each component, plus all rpointer.* files, into
  ``$DOUT_S_ROOT/rest/<date>/`` and then copies that whole directory back into
  the run directory. The run directory's restart files are therefore copies:
  a jump edits the copy, and the archived set stays pristine — it is the
  rollback source.
* Older restart files are deleted (or moved to ``<comp>/rest``) from the run
  directory at the next archive, so after the post-jump segment the run
  directory alone can no longer undo a jump when DOUT_S is on.

The SLURM probe mirrors exocam-casemgr's (job name == case name, one
``squeue --me`` snapshot); it is re-implemented here, not imported.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Set, Tuple

from .advise import ADVICE_SCHEMA_VERSION, TAPERED_ADVICE_SCHEMA_VERSION
from .ocean_advise import (COUPLED_SCHEMA_VERSION, OCEAN_ADVICE_SCHEMA_VERSION,
                           OCEAN_PATTERN_SCHEMA_VERSION, OCEAN_PROBE_SCHEMA_VERSION)
from .restart import (
    DATE_RE,
    STATE_DIR,
    backup_path,
    case_of,
    find_jump_logs,
    first_model_year,
    is_jumped,
    log_path,
    restart_date,
    restore,
)

_HIST_RE = re.compile(r"\.(\d{4})-(\d{2})(?:-\d{2}(?:-\d{5})?)?\.nc$")

#: Advice schema versions this build understands. ``preflight`` refuses advice
#: whose ``schema_version`` is missing or not in this set (feasibility-review
#: finding, Stage 0 item 5).
KNOWN_ADVICE_SCHEMA_VERSIONS = {ADVICE_SCHEMA_VERSION, TAPERED_ADVICE_SCHEMA_VERSION,
                                OCEAN_ADVICE_SCHEMA_VERSION, OCEAN_PATTERN_SCHEMA_VERSION,
                                OCEAN_PROBE_SCHEMA_VERSION, COUPLED_SCHEMA_VERSION}


def active_jobs() -> Optional[Set[str]]:
    """Names of the user's queued/running SLURM jobs; None if squeue is unusable."""
    try:
        res = subprocess.run(["squeue", "--me", "-h", "-o", "%j"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, timeout=60)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if res.returncode != 0:
        return None
    return {line.strip() for line in res.stdout.splitlines() if line.strip()}


JobProbe = Callable[[], Optional[Set[str]]]


@dataclass(frozen=True)
class Finding:
    level: str            # "ok" | "warn" | "block"
    message: str


def blocked(findings) -> bool:
    return any(f.level == "block" for f in findings)


def _job_finding(case: str, probe: Optional[JobProbe]) -> Finding:
    if probe is None:
        return Finding("warn", "SLURM check skipped (--skip-slurm-check): you have "
                               "confirmed the case is not queued or running")
    jobs = probe()
    if jobs is None:
        return Finding("block", "squeue unavailable: cannot confirm the case is not "
                                "queued or running (pass --skip-slurm-check only "
                                "after confirming it yourself)")
    if case in jobs:
        return Finding("block", f"a SLURM job named {case} is queued or running: "
                                f"stop the chain (or wait for the segment to end) "
                                f"before editing its restart files")
    return Finding("ok", f"no queued or running job named {case}")


def _ym(date: str) -> Tuple[int, int]:
    y, m, _, _ = DATE_RE.fullmatch(date).groups()
    return int(y), int(m)


def _hist_ym(name: str) -> Optional[Tuple[int, int]]:
    m = _HIST_RE.search(name)
    return (int(m.group(1)), int(m.group(2))) if m else None


def _rpointer_targets(rundir: Path) -> dict:
    """Map each ``rpointer.*`` file to the (first) restart file name it names."""
    out = {}
    for rp in sorted(rundir.glob("rpointer.*")):
        lines = [l.strip() for l in rp.read_text().splitlines()
                if l.strip() and not l.strip().startswith("#")]
        out[rp.name] = lines[0].lstrip("./").strip() if lines else None
    return out


def _rpointer_dates(rundir: Path) -> dict:
    out = {}
    for rp_name, name in _rpointer_targets(rundir).items():
        m = DATE_RE.search(name) if name else None
        out[rp_name] = "-".join(m.groups()) if m else None
    return out


def preflight(restart_file, archive: Optional[Path] = None,
              probe: Optional[JobProbe] = active_jobs,
              advice: Optional[dict] = None,
              allow_no_archive: bool = False) -> List[Finding]:
    """Checks before editing ``restart_file`` (a cice.r or docn.r) in place.
    Any "block" finding refuses.

    ``archive`` is the case's short-term archive root (DOUT_S_ROOT, i.e. the
    directory holding ``rest/``). ``probe`` = None skips the SLURM check.

    A verified rollback source is mandatory (feasibility-review finding 1):
    without ``archive``, ``preflight`` blocks unless ``allow_no_archive`` is
    also set, in which case every component restart named by a
    ``rpointer.*`` file must be present in the run directory itself (the
    rundir-sourced rollback path ``runstate.plan_rollback`` falls back to).
    """
    cice_r = Path(restart_file).resolve()
    rundir = cice_r.parent
    case = case_of(cice_r)
    date = restart_date(cice_r)
    out: List[Finding] = [_job_finding(case, probe)]

    # every component restarts from the same date as the file we edit
    dates = _rpointer_dates(rundir)
    wrong = {k: v for k, v in dates.items() if v != date}
    if not dates:
        out.append(Finding("block", f"no rpointer.* files in {rundir}"))
    elif wrong:
        out.append(Finding("block", f"rpointer files disagree with {cice_r.name} "
                                    f"({date}): {wrong} — stale or mixed pointers"))
    else:
        out.append(Finding("ok", f"all {len(dates)} rpointer files point at {date}"))

    # run directory must not hold output from beyond this restart
    h0 = [(_hist_ym(p.name), p.name) for p in rundir.glob(f"{case}.cam.h0.*.nc")]
    h0 = [x for x in h0 if x[0] is not None]
    if not h0:
        out.append(Finding("warn", "no cam.h0 files in the run directory (archived?): "
                                   "could not cross-check how far the run has got"))
    else:
        newest = max(h0)
        if newest[0] >= _ym(date):
            out.append(Finding("block", f"{newest[1]} is at or past the restart date "
                                        f"{date}: the run has moved on and the "
                                        f"rpointers are stale"))
        else:
            out.append(Finding("ok", f"newest history {newest[1]} precedes the restart"))

    # the archived set is the rollback source — a verified rollback source is
    # mandatory (feasibility-review finding 1): fail closed, do not warn.
    if archive is None:
        if not allow_no_archive:
            out.append(Finding("block", "no --archive given: a verified rollback "
                                        "source is required before a jump (pass "
                                        "--archive, or --allow-no-archive-rollback "
                                        "only after confirming every component "
                                        "restart is retained in the run directory)"))
        else:
            targets = _rpointer_targets(rundir)
            missing = {rp: name for rp, name in targets.items()
                      if name and not (rundir / name).is_file()}
            if not targets:
                out.append(Finding("block", f"--allow-no-archive-rollback: no "
                                            f"rpointer.* files in {rundir} to verify"))
            elif missing:
                out.append(Finding("block", f"--allow-no-archive-rollback: rpointer-"
                                            f"named restart(s) missing from {rundir}: "
                                            f"{missing} — rollback would be impossible"))
            else:
                out.append(Finding("ok", f"--allow-no-archive-rollback: all "
                                         f"{len(targets)} rpointer-named restarts "
                                         f"present in {rundir} (no archived copy "
                                         f"verified)"))
    else:
        rest = Path(archive) / "rest" / date
        if not (rest / cice_r.name).is_file():
            out.append(Finding("block", f"no archived restart set with {cice_r.name} "
                                        f"in {rest}: rollback after the next segment "
                                        f"would be impossible"))
        else:
            try:
                jumped = is_jumped(rest / cice_r.name)
            except ImportError:
                jumped = False
            if jumped:
                out.append(Finding("block", f"the archived {cice_r.name} is itself "
                                            f"jumped — it is not a pristine rollback copy"))
            else:
                out.append(Finding("ok", f"pristine archived restart set at {rest}"))

    # the advice must be bound to this case, this restart date, and a schema
    # this build understands (feasibility-review finding, Stage 0 item 5)
    if advice is not None:
        schema = advice.get("schema_version")
        if schema not in KNOWN_ADVICE_SCHEMA_VERSIONS:
            out.append(Finding("block", f"advice schema_version {schema!r} is missing "
                                        f"or unrecognized (known: "
                                        f"{sorted(KNOWN_ADVICE_SCHEMA_VERSIONS)}): "
                                        f"refusing to bind it to this jump"))
        if not advice.get("case"):
            out.append(Finding("block", "advice has no case name recorded: cannot "
                                        "verify it is bound to this restart"))
        elif advice["case"] != case:
            out.append(Finding("block", f"advice is for {advice['case']!r}, the "
                                        f"restart is {case!r}"))
        yr = advice.get("model_year")
        if yr is None:
            out.append(Finding("block", "advice has no model_year recorded: cannot "
                                        "verify it is bound to this restart date"))
        else:
            gap = first_model_year(date) - 1 - int(yr)
            if gap > 5:
                out.append(Finding("block", f"advice uses data through model year "
                                            f"{yr} but the restart is {date}: "
                                            f"{gap} years stale — re-run advise"))
            elif gap > 0:
                out.append(Finding("warn", f"advice is {gap} year(s) older than the "
                                           f"restart; fine for the factor, the landing "
                                           f"check allows for it"))
            elif gap < 0:
                out.append(Finding("block", f"advice uses data through model year {yr}, "
                                            f"after this restart ({date})"))
    return out


# ---------------------------------------------------------------------------
# rollback
# ---------------------------------------------------------------------------

@dataclass
class RollbackPlan:
    rundir: Path
    case: str
    date: str
    source: str                          # "archive" | "rundir"
    findings: List[Finding] = field(default_factory=list)
    actions: List[str] = field(default_factory=list)
    discarded: List[str] = field(default_factory=list)
    executed: bool = False


#: which file's log speaks for a jump that edited several files at one date
#: (a coupled ocean + atmosphere jump): the one ``check`` scores
_PRIMARY = {"som_ocean": 0, "aqua_ice": 1, "cam_atm": 2}


def jump_log_group(rundir: Path, date: Optional[str]) -> List[Path]:
    """Active jump logs of one jump (all files edited at one restart date),
    primary first. Several dates active without ``date`` is an error."""
    logs = find_jump_logs(rundir)
    if date:
        logs = [p for p in logs if date in p.name]
    if not logs:
        raise FileNotFoundError(f"no active jump log in {rundir / STATE_DIR}"
                                + (f" for {date}" if date else ""))
    metas = {p: json.loads(p.read_text()) for p in logs}
    dates = sorted({m.get("restart_date") for m in metas.values()})
    if len(dates) > 1:
        raise ValueError(f"active jump logs at several dates {dates}, pass --date")
    return sorted(logs, key=lambda p: _PRIMARY.get(metas[p].get("plugin"), 9))


def select_jump_log(rundir: Path, date: Optional[str]) -> Path:
    return jump_log_group(rundir, date)[0]


def _after(name: str, date: str) -> bool:
    """Output produced by the discarded segment: restarts dated after the jump,
    history for months at or after it."""
    m = DATE_RE.findall(name)
    if m:
        return "-".join(m[-1]) > date
    ym = _hist_ym(name)
    return ym is not None and ym >= _ym(date)


def plan_rollback(rundir, archive: Optional[Path] = None, date: Optional[str] = None,
                  probe: Optional[JobProbe] = active_jobs) -> RollbackPlan:
    """Plan resetting a case to its pre-jump restart set (nothing is changed)."""
    rundir = Path(rundir).resolve()
    log = select_jump_log(rundir, date)
    meta = json.loads(log.read_text())
    date, case = meta["restart_date"], meta["case"]
    cice = rundir / meta["restart_file"]
    rest = Path(archive) / "rest" / date if archive is not None else None
    use_archive = rest is not None and rest.is_dir()
    plan = RollbackPlan(rundir, case, date, "archive" if use_archive else "rundir")
    plan.findings.append(_job_finding(case, probe))

    if use_archive:
        files = sorted(p for p in rest.iterdir() if p.is_file())
        arch_cice = rest / cice.name
        if not arch_cice.is_file():
            plan.findings.append(Finding("block", f"{rest} has no {cice.name}"))
        else:
            try:
                if is_jumped(arch_cice):
                    plan.findings.append(Finding("block", f"archived {cice.name} is "
                                                          f"jumped, not pristine"))
            except ImportError:
                plan.findings.append(Finding("warn", "netCDF4 unavailable: archived "
                                                     "cice.r not verified pristine"))
        plan.actions.append(f"save current rpointer.* to {STATE_DIR}/rollback-{date}/")
        plan.actions += [f"copy {rest.name}/{p.name} -> run/{p.name}" for p in files]
    else:
        if archive is not None:
            plan.findings.append(Finding("warn", f"no archived set {rest}; rolling "
                                                 f"back from the run directory"))
        if not backup_path(cice).exists():
            plan.findings.append(Finding("block", f"no pristine backup "
                                                  f"{backup_path(cice)} and no archived "
                                                  f"set: cannot roll back cice.r"))
        plan.actions.append(f"save current rpointer.* to {STATE_DIR}/rollback-{date}/")
        for rp in sorted(rundir.glob("rpointer.*")):
            text = rp.read_text()
            new = DATE_RE.sub(date, text)
            for line in new.splitlines():
                name = line.strip().lstrip("./").strip()
                if not name or name.startswith("#"):
                    continue
                if not (rundir / name).is_file():
                    plan.findings.append(Finding("block", f"{rp.name} would name "
                                                          f"{name}, missing from the "
                                                          f"run directory (archived "
                                                          f"away? pass --archive)"))
            if new != text:
                plan.actions.append(f"rewrite {rp.name} to {date}")
        for lg in jump_log_group(rundir, date):
            f = rundir / json.loads(lg.read_text())["restart_file"]
            if not backup_path(f).exists():
                plan.findings.append(Finding("block", f"no pristine backup "
                                                      f"{backup_path(f)} for {f.name}"))
            plan.actions.append(f"restore {f.name} from its pristine backup")
    for lg in jump_log_group(rundir, date):
        plan.actions.append(f"retire jump log {lg.name} (-> .rolledback.json)")

    for p in sorted(rundir.iterdir()):
        if p.is_file() and p.name.startswith(case + ".") and _after(p.name, date):
            plan.discarded.append(str(p))
    if archive is not None:
        for p in sorted(Path(archive).glob("*/hist/*")):
            if p.name.startswith(case + ".") and _after(p.name, date):
                plan.discarded.append(str(p))
    return plan


def execute_rollback(plan: RollbackPlan, archive: Optional[Path] = None) -> RollbackPlan:
    if blocked(plan.findings):
        raise RuntimeError("rollback plan has blocking findings")
    rundir = plan.rundir
    save = rundir / STATE_DIR / f"rollback-{plan.date}"
    save.mkdir(parents=True, exist_ok=True)
    for rp in rundir.glob("rpointer.*"):
        shutil.copy2(rp, save / rp.name)

    group = jump_log_group(rundir, plan.date)
    if plan.source == "archive":
        rest = Path(archive) / "rest" / plan.date
        for p in rest.iterdir():
            if p.is_file():
                shutil.copy2(p, rundir / p.name)
        for log in group:
            log.rename(log.with_name(log.name[: -len(".json")] + ".rolledback.json"))
    else:
        for rp in rundir.glob("rpointer.*"):
            rp.write_text(DATE_RE.sub(plan.date, rp.read_text()))
        for log in group:
            restore(rundir / json.loads(log.read_text())["restart_file"],
                    retire_log="rolledback")
    plan.executed = True
    return plan
