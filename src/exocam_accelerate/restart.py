"""Restart-file layer: apply an aqua_ice jump to cice.r in place, with backup.

The only module that imports netCDF4 (optional extra ``[netcdf]``).

Mechanism (decided 2026-09-25, docs/restart-integration-questions.md §7):
in-place edit + continuation. The run directory's rpointer.ice keeps naming
the same file, so the user simply resubmits with CONTINUE_RUN=TRUE.

Backup discipline:
* The first jump on a file copies it to ``<stem>.pre-accel.nc`` (pristine).
* Any later jump on the same file re-reads from the pristine copy, so jumps
  never compound by re-running the command; ``restore`` puts it back.
* The jumped file carries a global attribute ``exocam_accelerate`` (JSON)
  and a sidecar ``<file>.accel.json`` log; a file that carries the attribute
  but has no pristine backup is refused rather than scaled twice.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .aqua_ice import ICE_FIELDS, SNOW_FIELDS, AquaIcePlugin

try:
    import netCDF4
except ImportError:  # pragma: no cover - exercised only without the extra
    netCDF4 = None

ATTR = "exocam_accelerate"
BACKUP_SUFFIX = ".pre-accel.nc"
LOG_SUFFIX = ".accel.json"


def _require_netcdf():
    if netCDF4 is None:
        raise ImportError("the restart layer needs netCDF4: "
                          "pip install 'exocam-accelerate[netcdf]'")


def backup_path(path) -> Path:
    path = Path(path)
    if path.name.endswith(BACKUP_SUFFIX):
        raise ValueError(f"{path} is itself a backup")
    return path.with_name(path.name[:-3] + BACKUP_SUFFIX if path.suffix == ".nc"
                          else path.name + BACKUP_SUFFIX)


def log_path(path) -> Path:
    path = Path(path)
    return path.with_name(path.name + LOG_SUFFIX)


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


def _has_attr(path) -> bool:
    with netCDF4.Dataset(path, "r") as ds:
        return ATTR in ds.ncattrs()


@dataclass(frozen=True)
class JumpRecord:
    path: Path
    backup: Path
    ice_factor: float
    snow_factor: float
    written: bool
    adjustments: Dict[str, str]
    metadata: dict


def apply_ice_jump(path, ice_factor: float, snow_factor: float = 1.0,
                   provenance: Optional[dict] = None,
                   dry_run: bool = False) -> JumpRecord:
    """Scale vicen/eicen (and optionally vsnon/esnon) in a cice.r file.

    ``provenance`` (e.g. ``Advice.to_dict()``) is stored in the sidecar log so
    the jump can be audited and the post-jump run checked against it.
    """
    _require_netcdf()
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    backup = backup_path(path)

    if backup.exists():
        if _has_attr(backup):
            raise RuntimeError(f"{backup} is marked as already jumped; it is not "
                               f"a pristine backup — refusing")
        source = backup
    else:
        if _has_attr(path):
            raise RuntimeError(f"{path} was already jumped but its pristine backup "
                               f"{backup.name} is missing — refusing to scale twice")
        source = path

    plugin = AquaIcePlugin()
    names = ICE_FIELDS + (SNOW_FIELDS if snow_factor != 1.0 else ())
    before = _read(source, names)
    jumped = plugin.apply_delta(before, (ice_factor, snow_factor))
    after, report = plugin.enforce_constraints(before, jumped)

    meta = {
        "tool": "exocam-accelerate",
        "plugin": plugin.name,
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "ice_factor": float(ice_factor),
        "snow_factor": float(snow_factor),
        "fields": list(after),
        "pristine_backup": backup.name,
    }
    if dry_run:
        return JumpRecord(path, backup, ice_factor, snow_factor, False,
                          dict(report.adjustments), meta)

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

    log = dict(meta, adjustments=dict(report.adjustments), advice=provenance)
    log_path(path).write_text(json.dumps(log, indent=2))
    return JumpRecord(path, backup, ice_factor, snow_factor, True,
                      dict(report.adjustments), meta)


def restore(path) -> Path:
    """Put the pristine copy back over a jumped file (backup is kept)."""
    path = Path(path).resolve()
    backup = backup_path(path)
    if not backup.exists():
        raise FileNotFoundError(f"no pristine backup {backup.name} for {path.name}")
    shutil.copy2(backup, path)
    log = log_path(path)
    if log.exists():
        log.unlink()
    return backup
