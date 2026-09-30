"""Restart-file layer: apply an aqua_ice jump to cice.r in place, with backup.

The only module that imports netCDF4 (optional extra ``[netcdf]``).

Mechanism (decided 2026-09-25, docs/restart-integration-questions.md §7):
in-place edit + continuation. rpointer.ice keeps naming the same file, so the
case is simply resubmitted with CONTINUE_RUN=TRUE.

Where the bookkeeping lives — ``<rundir>/exocam_accelerate/``:
* ``<cice.r name>.pre-accel.nc`` — pristine copy made on the first jump;
  re-running a jump re-reads it, so jumps never compound.
* ``<cice.r name>.accel.json`` — the jump log (factor, restart date, first
  post-jump model year, full advice); the input to ``check`` and ``rollback``.
They must NOT sit next to the restart: CESM 1.2's st_archive.sh globs
``${CASE}.cice.r.[0-9]*`` in the run directory at the end of every segment and
deletes (or moves) everything but the newest match — a backup or log named
after the restart would be swept away before it is needed. st_archive never
descends into subdirectories.

A tapered jump (``taper.py``) is sized by ``taper`` from the case's archived
restarts and a companion ``<advice>.taper.nc`` weight map; ``jump`` applies
``1 + (F-1)*weight`` per cell and keeps a copy of the map as
``<cice.r name>.taper.nc`` beside the log.

The jumped file also carries a JSON global attribute ``exocam_accelerate``; a
file that has it but no pristine backup is refused rather than scaled twice.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .aqua_ice import ICE_FIELDS, SNOW_FIELDS, AquaIcePlugin, check_factor_map
from .taper import TaperConfig, TaperMask, cell_factors, stefan_mask

try:
    import netCDF4
except ImportError:  # pragma: no cover - exercised only without the extra
    netCDF4 = None

ATTR = "exocam_accelerate"
STATE_DIR = "exocam_accelerate"
BACKUP_SUFFIX = ".pre-accel.nc"
LOG_SUFFIX = ".accel.json"
TAPER_SUFFIX = ".taper.nc"
DATE_RE = re.compile(r"(\d{4})-(\d{2})-(\d{2})-(\d{5})")


def _require_netcdf():
    if netCDF4 is None:
        raise ImportError("the restart layer needs netCDF4: "
                          "pip install 'exocam-accelerate[netcdf]'")


def state_dir(path) -> Path:
    """Bookkeeping directory for a restart file in a run directory."""
    return Path(path).resolve().parent / STATE_DIR


def backup_path(path) -> Path:
    name = Path(path).name
    stem = name[:-3] if name.endswith(".nc") else name
    return state_dir(path) / (stem + BACKUP_SUFFIX)


def log_path(path) -> Path:
    return state_dir(path) / (Path(path).name + LOG_SUFFIX)


def restart_date(path) -> str:
    """``YYYY-MM-DD-SSSSS`` from a CESM restart file name."""
    m = DATE_RE.findall(Path(path).name)
    if not m:
        raise ValueError(f"no restart date in {Path(path).name}")
    return "-".join(m[-1])


def case_of(path) -> str:
    name = Path(path).name
    if ".cice.r." not in name:
        raise ValueError(f"{name} is not a cice.r file")
    return name.split(".cice.r.")[0]


def first_model_year(date: str) -> int:
    """First complete model year run from a restart written at ``date``."""
    y, m, d, s = DATE_RE.fullmatch(date).groups()
    return int(y) if (m, d, s) == ("01", "01", "00000") else int(y) + 1


def locate_cice_restart(rundir) -> Path:
    """The cice.r file rpointer.ice names — what a continuation will read."""
    rundir = Path(rundir)
    rp = rundir / "rpointer.ice"
    if not rp.is_file():
        raise FileNotFoundError(f"no rpointer.ice in {rundir}")
    name = rp.read_text().split()[0]
    path = (rundir / name).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"rpointer.ice names {name}, which does not exist")
    if ".cice.r." not in path.name:
        raise ValueError(f"rpointer.ice names {path.name}, not a cice.r file")
    return path


def _read(path, names) -> Dict[str, np.ndarray]:
    with netCDF4.Dataset(path, "r") as ds:
        ds.set_auto_mask(False)
        missing = [n for n in names if n not in ds.variables]
        if missing:
            raise KeyError(f"{path}: missing variables {missing}")
        return {n: np.array(ds.variables[n][:]) for n in names}


def is_jumped(path) -> bool:
    """True when the file carries the jump attribute (header read only)."""
    _require_netcdf()
    with netCDF4.Dataset(path, "r") as ds:
        return ATTR in ds.ncattrs()


@dataclass(frozen=True)
class JumpRecord:
    path: Path
    backup: Path
    log: Path
    ice_factor: float
    snow_factor: float
    written: bool
    adjustments: Dict[str, str]
    metadata: dict


def apply_ice_jump(path, ice_factor: float, snow_factor: float = 1.0,
                   provenance: Optional[dict] = None,
                   dry_run: bool = False,
                   taper_file=None) -> JumpRecord:
    """Scale vicen/eicen (and optionally vsnon/esnon) in a cice.r file.

    ``provenance`` (e.g. ``Advice.to_dict()``) is stored in the jump log so the
    post-jump run can be checked against it.

    With ``taper_file`` (from ``write_taper_file``), ``ice_factor`` is the peak
    factor and each cell gets ``1 + (ice_factor - 1) * weight``; the map is
    copied next to the log.
    """
    _require_netcdf()
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    backup, log = backup_path(path), log_path(path)

    if backup.exists():
        if is_jumped(backup):
            raise RuntimeError(f"{backup} is marked as already jumped; it is not "
                               f"a pristine backup — refusing")
        source = backup
    else:
        if is_jumped(path):
            raise RuntimeError(f"{path.name} was already jumped but its pristine "
                               f"backup {backup} is missing — refusing to scale twice")
        source = path

    plugin = AquaIcePlugin()
    names = ICE_FIELDS + (SNOW_FIELDS if snow_factor != 1.0 else ())
    before = _read(source, names)
    ice_f = ice_factor
    taper_copy = None
    if taper_file is not None:
        weight = read_taper_file(taper_file)["weight"]
        if weight.shape != before["vicen"].shape[1:]:
            raise ValueError(f"taper map {Path(taper_file).name} is {weight.shape}, "
                             f"the restart grid is {before['vicen'].shape[1:]}")
        ice_f = check_factor_map(cell_factors(ice_factor, weight), "ice factor")
        taper_copy = state_dir(path) / (path.name + TAPER_SUFFIX)
    jumped = plugin.apply_delta(before, (ice_f, snow_factor))
    after, report = plugin.enforce_constraints(before, jumped)

    date = restart_date(path)
    meta = {
        "tool": "exocam-accelerate",
        "plugin": plugin.name,
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "case": case_of(path),
        "restart_file": path.name,
        "restart_date": date,
        "jump_model_year": first_model_year(date),
        "ice_factor": float(ice_factor),
        "snow_factor": float(snow_factor),
        "fields": list(after),
        "pristine_backup": str(backup),
    }
    if taper_file is not None:
        meta["taper"] = {"source": str(Path(taper_file).resolve()),
                         "sha256": file_sha256(taper_file),
                         "copy": str(taper_copy),
                         "cell_factor_min": float(np.min(ice_f)),
                         "cell_factor_max": float(np.max(ice_f))}
    if dry_run:
        return JumpRecord(path, backup, log, ice_factor, snow_factor, False,
                          dict(report.adjustments), meta)

    backup.parent.mkdir(exist_ok=True)
    if taper_copy is not None:
        shutil.copy2(taper_file, taper_copy)
    if source is path:
        shutil.copy2(path, backup)
    else:
        shutil.copy2(backup, path)

    with netCDF4.Dataset(path, "r+") as ds:
        ds.set_auto_mask(False)
        for name, arr in after.items():
            var = ds.variables[name]
            var[:] = arr.astype(var.dtype)
        ds.setncattr(ATTR, json.dumps(meta))

    check = _read(path, list(after))
    for name, arr in after.items():
        if not np.array_equal(check[name], arr.astype(check[name].dtype)):
            raise RuntimeError(f"verification failed for {name} in {path}; "
                               f"restore with the pristine copy {backup}")

    log.write_text(json.dumps(dict(meta, adjustments=dict(report.adjustments),
                                   advice=provenance), indent=2))
    return JumpRecord(path, backup, log, ice_factor, snow_factor, True,
                      dict(report.adjustments), meta)


def restore(path, retire_log: str = "restored") -> Path:
    """Put the pristine copy back over a jumped file (backup is kept).

    The jump log is kept, renamed ``*.accel.<retire_log>.json``, so the history
    of what was tried survives.
    """
    path = Path(path).resolve()
    backup = backup_path(path)
    if not backup.exists():
        raise FileNotFoundError(f"no pristine backup {backup} for {path.name}")
    shutil.copy2(backup, path)
    log = log_path(path)
    if log.exists():
        log.rename(log.with_name(log.name[: -len(".json")] + f".{retire_log}.json"))
    return backup


def find_jump_logs(rundir) -> list:
    """Active (not restored / rolled back) jump logs in a run directory."""
    d = Path(rundir) / STATE_DIR
    return sorted(p for p in d.glob(f"*{LOG_SUFFIX}")
                  if p.name.endswith(".nc" + LOG_SUFFIX))


# ---------------------------------------------------------------------------
# tapered jump inputs (taper.py)
# ---------------------------------------------------------------------------

def file_sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def shift_date(date: str, years: int) -> str:
    y, m, d, s = DATE_RE.fullmatch(date).groups()
    if int(y) - years < 1:
        raise ValueError(f"{date} minus {years} years is before model year 1")
    return f"{int(y) - years:04d}-{m}-{d}-{s}"


def archived_cice_r(archive, case: str, date: str) -> Path:
    p = Path(archive) / "rest" / date / f"{case}.cice.r.{date}.nc"
    if not p.is_file():
        raise FileNotFoundError(f"no archived cice.r for {date}: {p}")
    return p


def find_grid_file(archive, case: str) -> Path:
    """Any cice.h history file of the case — the grid (tarea, tmask) source."""
    d = Path(archive) / "ice" / "hist"
    files = sorted(d.glob(f"{case}.cice.h.*.nc"))
    if not files:
        raise FileNotFoundError(f"no {case}.cice.h.*.nc in {d} (pass --grid-file)")
    return files[-1]


def read_ice_state(path):
    """Grid-box ice volume per area (sum of vicen) and ice fraction (sum of aicen)."""
    _require_netcdf()
    f = _read(path, ["vicen", "aicen"])
    return f["vicen"].sum(axis=0), f["aicen"].sum(axis=0)


def read_cell_area(grid_file) -> np.ndarray:
    """tarea masked by tmask from a cice.h file."""
    _require_netcdf()
    g = _read(grid_file, ["tarea", "tmask"])
    return np.where(g["tmask"] > 0, g["tarea"], 0.0).astype(float)


def build_taper_mask(archive, case: str, date: str, baseline_years: int = 10,
                     grid_file=None, config: TaperConfig = TaperConfig()):
    """Stefan mask from the pristine archived restarts at ``date`` and
    ``baseline_years`` earlier. Returns ``(mask, sources)``."""
    if baseline_years <= 0:
        raise ValueError(f"baseline_years must be positive, got {baseline_years!r}")
    now = archived_cice_r(archive, case, date)
    old = archived_cice_r(archive, case, shift_date(date, baseline_years))
    for p in (now, old):
        if is_jumped(p):
            raise RuntimeError(f"{p.name} is marked as jumped: the taper baseline "
                               f"must be pristine archived restarts")
    grid = Path(grid_file) if grid_file else find_grid_file(archive, case)
    h_now, aice = read_ice_state(now)
    h_old, _ = read_ice_state(old)
    area = read_cell_area(grid)
    if area.shape != h_now.shape:
        raise ValueError(f"grid file {grid.name} is {area.shape}, restarts are "
                         f"{h_now.shape}")
    mask = stefan_mask(h_old, h_now, float(baseline_years), aice, area, config)
    sources = {str(p): file_sha256(p) for p in (old, now, grid)}
    return mask, sources


TAPER_FIELDS = ("weight", "w", "growth", "h_now", "area")


def write_taper_file(path, mask: TaperMask, peak: Optional[float] = None) -> str:
    """Save the mask (and, with ``peak``, the cell factors); returns its sha256."""
    _require_netcdf()
    with netCDF4.Dataset(path, "w") as ds:
        nj, ni = mask.weight.shape
        ds.createDimension("nj", nj)
        ds.createDimension("ni", ni)
        units = {"weight": "1", "w": "1", "growth": "m/yr", "h_now": "m",
                 "area": "m2"}
        for name in TAPER_FIELDS:
            v = ds.createVariable(name, "f8", ("nj", "ni"))
            v[:] = getattr(mask, name)
            v.units = units[name]
        if peak is not None:
            v = ds.createVariable("ice_factor", "f8", ("nj", "ni"))
            v[:] = cell_factors(peak, mask.weight)
            v.units = "1"
        ds.setncattr(ATTR, json.dumps(dict(mask.summary(), peak_factor=peak)))
    return file_sha256(path)


def read_taper_file(path) -> Dict[str, np.ndarray]:
    _require_netcdf()
    return _read(path, list(TAPER_FIELDS))
