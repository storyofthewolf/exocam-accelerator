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


#: restart kinds this layer edits in place
KINDS = (".cice.r.", ".docn.r.", ".cam.r.")


def case_of(path) -> str:
    name = Path(path).name
    for kind in KINDS:
        if kind in name:
            return name.split(kind)[0]
    raise ValueError(f"{name} is not a cice.r, docn.r or cam.r file")


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


def locate_docn_restart(rundir) -> Path:
    """The docn.r file rpointer.ocn names (its first line; the second is the
    stream restart docn.rs1.bin)."""
    rundir = Path(rundir)
    rp = rundir / "rpointer.ocn"
    if not rp.is_file():
        raise FileNotFoundError(f"no rpointer.ocn in {rundir}")
    lines = [l.strip() for l in rp.read_text().splitlines() if l.strip()]
    if not lines:
        raise ValueError(f"{rp} is empty")
    path = (rundir / lines[0]).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"rpointer.ocn names {lines[0]}, which does not exist")
    if ".docn.r." not in path.name:
        raise ValueError(f"rpointer.ocn names {path.name}, not a docn.r file")
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


# ---------------------------------------------------------------------------
# som_ocean: docn.r somtp (som_ocean.py)
# ---------------------------------------------------------------------------

PATTERN_SUFFIX = ".pattern.nc"


def read_somtp(path) -> np.ndarray:
    _require_netcdf()
    return _read(path, ["somtp"])["somtp"].astype(float)


def find_docn_domain(rundir) -> Optional[Path]:
    """The docn domain file named in the run directory's ``docn_ocn_in``
    (``domainfile = '...'``); None when absent or unreadable."""
    nml = Path(rundir) / "docn_ocn_in"
    if not nml.is_file():
        return None
    m = re.search(r"domainfile\s*=\s*['\"]([^'\"]+)['\"]", nml.read_text())
    if not m:
        return None
    p = Path(m.group(1))
    return p if p.is_file() else None


def read_ocean_grid(domain_file) -> Dict[str, np.ndarray]:
    """``mask``, ``area``, ``lat``, ``lon`` as (nj, ni) from a docn domain file.

    Accepts a CESM domain file (``xc``/``yc``/``mask``/``area`` on (nj, ni)) or
    a SOM forcing file like ``pop_frc`` (1-D ``xc(lon)``, ``yc(lat)``). The
    grid must be the docn domain: its nj*ni must equal docn.r's gsize.
    """
    _require_netcdf()
    with netCDF4.Dataset(domain_file) as ds:
        ds.set_auto_mask(False)
        missing = [n for n in ("mask", "area", "xc", "yc") if n not in ds.variables]
        if missing:
            raise KeyError(f"{domain_file}: missing {missing} (need a docn domain "
                           f"or pop_frc file)")
        mask = np.array(ds.variables["mask"][:], dtype=float)
        area = np.array(ds.variables["area"][:], dtype=float)
        xc = np.array(ds.variables["xc"][:], dtype=float)
        yc = np.array(ds.variables["yc"][:], dtype=float)
    mask, area = np.squeeze(mask), np.squeeze(area)
    if mask.ndim != 2:
        raise ValueError(f"{domain_file}: mask is {mask.shape}, expected 2-D")
    nj, ni = mask.shape
    xc, yc = np.squeeze(xc), np.squeeze(yc)
    if xc.ndim == 1 and yc.ndim == 1:
        yc, xc = np.meshgrid(yc, xc, indexing="ij")
    if xc.shape != (nj, ni) or yc.shape != (nj, ni):
        raise ValueError(f"{domain_file}: coordinates {xc.shape}/{yc.shape} do not "
                         f"match the {nj}x{ni} mask")
    return {"mask": mask, "area": area, "lat": yc, "lon": xc}


def apply_ocean_jump(path, somtp_dT: float, provenance: Optional[dict] = None,
                     dry_run: bool = False, pattern_file=None,
                     domain_file=None) -> JumpRecord:
    """Add ``somtp_dT`` (K) to somtp in a docn.r file, in place, with backup.

    With ``pattern_file`` (from ``write_pattern_file``) each cell gets
    ``somtp_dT * weight`` (area mean ``somtp_dT``); the map is copied beside
    the log. ``domain_file`` supplies the ocean mask (masked cells untouched)
    and the area weights for the report; without it every non-frozen cell
    gets the increment. Same backup / log / re-apply-from-pristine rules as
    ``apply_ice_jump``; the record's ``ice_factor`` slot holds the increment.
    """
    from .som_ocean import SomOceanPlugin, TK_FRZ_SW, area_mean, check_dT

    _require_netcdf()
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if ".docn.r." not in path.name:
        raise ValueError(f"{path.name} is not a docn.r file")
    backup, log = backup_path(path), log_path(path)
    if backup.exists():
        if is_jumped(backup):
            raise RuntimeError(f"{backup} is marked as already jumped; it is not "
                               f"a pristine backup — refusing")
        source = backup
    else:
        if is_jumped(path):
            raise RuntimeError(f"{path.name} was already jumped but its pristine "
                               f"backup {backup} is missing — refusing to shift twice")
        source = path

    check_dT(somtp_dT)
    before = {"somtp": read_somtp(source)}
    n = before["somtp"].size
    ocean = area = None
    if domain_file is not None:
        g = read_ocean_grid(domain_file)
        if g["mask"].size != n:
            raise ValueError(f"domain {Path(domain_file).name} has {g['mask'].size} "
                             f"cells, somtp has {n}")
        ocean = g["mask"].ravel() > 0
        area = g["area"].ravel()
    delta = float(somtp_dT)
    pattern_copy = None
    if pattern_file is not None:
        weight = read_pattern_file(pattern_file)["weight"].ravel()
        if weight.size != n:
            raise ValueError(f"pattern map {Path(pattern_file).name} has "
                             f"{weight.size} cells, somtp has {n}")
        delta = check_dT(somtp_dT * weight, "per-cell somtp increment")
        pattern_copy = state_dir(path) / (path.name + PATTERN_SUFFIX)

    plugin = SomOceanPlugin()
    jumped = plugin.apply_delta(before, delta)
    after, report = plugin.enforce_constraints(before, jumped, ocean)
    d = after["somtp"] - before["somtp"]
    open_sea = before["somtp"] > TK_FRZ_SW + 0.05
    if ocean is not None:
        open_sea &= ocean
    mean_dT = (area_mean(d, area, open_sea) if area is not None
               else float(d[open_sea].mean()) if open_sea.any() else 0.0)

    date = restart_date(path)
    meta = {
        "tool": "exocam-accelerate",
        "plugin": plugin.name,
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "case": case_of(path),
        "restart_file": path.name,
        "restart_date": date,
        "jump_model_year": first_model_year(date),
        "somtp_dT": float(somtp_dT),
        "somtp_dT_applied_mean": float(mean_dT),
        "cell_dT_min": float(d.min()),
        "cell_dT_max": float(d.max()),
        "fields": ["somtp"],
        "pristine_backup": str(backup),
        "domain_file": str(Path(domain_file).resolve()) if domain_file else None,
    }
    if pattern_file is not None:
        meta["pattern"] = {"source": str(Path(pattern_file).resolve()),
                           "sha256": file_sha256(pattern_file),
                           "copy": str(pattern_copy)}
    if dry_run:
        return JumpRecord(path, backup, log, float(somtp_dT), 1.0, False,
                          dict(report.adjustments), meta)

    backup.parent.mkdir(exist_ok=True)
    if pattern_copy is not None:
        shutil.copy2(pattern_file, pattern_copy)
    if source is path:
        shutil.copy2(path, backup)
    else:
        shutil.copy2(backup, path)
    with netCDF4.Dataset(path, "r+") as ds:
        ds.set_auto_mask(False)
        var = ds.variables["somtp"]
        var[:] = after["somtp"].astype(var.dtype)
        ds.setncattr(ATTR, json.dumps(meta))
    check = read_somtp(path)
    if not np.array_equal(check, after["somtp"].astype(check.dtype)):
        raise RuntimeError(f"verification failed for somtp in {path}; restore with "
                           f"the pristine copy {backup}")
    log.write_text(json.dumps(dict(meta, adjustments=dict(report.adjustments),
                                   advice=provenance), indent=2))
    return JumpRecord(path, backup, log, float(somtp_dT), 1.0, True,
                      dict(report.adjustments), meta)


def archived_docn_r(archive, case: str, date: str) -> Path:
    p = Path(archive) / "rest" / date / f"{case}.docn.r.{date}.nc"
    if not p.is_file():
        raise FileNotFoundError(f"no archived docn.r for {date}: {p}")
    return p


def build_ocean_pattern(archive, case: str, date: str, domain_file,
                        baseline_years: int = 10, config=None):
    """Warming pattern from the pristine archived docn.r at ``date`` and
    ``baseline_years`` earlier. Returns ``(pattern, grid, sources)``."""
    from .som_ocean import PatternConfig, pattern_weights, to_grid

    if baseline_years <= 0:
        raise ValueError(f"baseline_years must be positive, got {baseline_years!r}")
    now = archived_docn_r(archive, case, date)
    old = archived_docn_r(archive, case, shift_date(date, baseline_years))
    for p in (now, old):
        if is_jumped(p):
            raise RuntimeError(f"{p.name} is marked as jumped: the pattern baseline "
                               f"must be pristine archived restarts")
    grid = read_ocean_grid(domain_file)
    nj, ni = grid["mask"].shape
    t_now = to_grid(read_somtp(now), nj, ni)
    t_old = to_grid(read_somtp(old), nj, ni)
    area = np.where(grid["mask"] > 0, grid["area"], 0.0)
    pat = pattern_weights(t_old, t_now, float(baseline_years), area,
                          grid["mask"] > 0, config or PatternConfig())
    sources = {str(p): file_sha256(p) for p in (old, now, Path(domain_file))}
    return pat, grid, sources


PATTERN_FIELDS = ("weight", "rate", "somtp_now", "area", "lat", "lon")


def write_pattern_file(path, pattern, grid, somtp_dT: Optional[float] = None) -> str:
    """Save the pattern on (nj, ni) with lat/lon (and, with ``somtp_dT``, the
    per-cell increment); returns its sha256."""
    _require_netcdf()
    with netCDF4.Dataset(path, "w") as ds:
        nj, ni = pattern.weight.shape
        ds.createDimension("nj", nj)
        ds.createDimension("ni", ni)
        vals = {"weight": (pattern.weight, "1"), "rate": (pattern.rate, "K/yr"),
                "somtp_now": (pattern.somtp_now, "K"), "area": (pattern.area, "1"),
                "lat": (grid["lat"], "degrees_north"), "lon": (grid["lon"], "degrees_east")}
        for name, (arr, units) in vals.items():
            v = ds.createVariable(name, "f8", ("nj", "ni"))
            v[:] = arr
            v.units = units
        if somtp_dT is not None:
            v = ds.createVariable("somtp_dT", "f8", ("nj", "ni"))
            v[:] = somtp_dT * pattern.weight
            v.units = "K"
        ds.setncattr(ATTR, json.dumps(dict(pattern.summary(), somtp_dT=somtp_dT)))
    return file_sha256(path)


def read_pattern_file(path) -> Dict[str, np.ndarray]:
    _require_netcdf()
    return _read(path, list(PATTERN_FIELDS))


def write_somtp_map(docn_r, domain_file, out) -> Path:
    """docn.r somtp -> a lat-lon netCDF (``somtp(lat, lon)`` in K, plus
    ``mask``/``area``) that ncview, Panoply or the viewer can show. For a
    rectilinear docn grid the coordinates are written 1-D."""
    from .som_ocean import to_grid

    grid = read_ocean_grid(domain_file)
    nj, ni = grid["mask"].shape
    somtp = to_grid(read_somtp(docn_r), nj, ni)
    rect = (np.allclose(grid["lat"], grid["lat"][:, :1]) and
            np.allclose(grid["lon"], grid["lon"][:1, :]))
    with netCDF4.Dataset(out, "w") as ds:
        if rect:
            ds.createDimension("lat", nj)
            ds.createDimension("lon", ni)
            la = ds.createVariable("lat", "f8", ("lat",))
            la[:] = grid["lat"][:, 0]
            la.units = "degrees_north"
            lo = ds.createVariable("lon", "f8", ("lon",))
            lo[:] = grid["lon"][0, :]
            lo.units = "degrees_east"
            dims = ("lat", "lon")
        else:
            ds.createDimension("nj", nj)
            ds.createDimension("ni", ni)
            dims = ("nj", "ni")
            for name, units in (("lat", "degrees_north"), ("lon", "degrees_east")):
                v = ds.createVariable(name, "f8", dims)
                v[:] = grid[name]
                v.units = units
        for name, arr, units in (("somtp", somtp, "K"), ("mask", grid["mask"], "1"),
                                 ("area", grid["area"], "1")):
            v = ds.createVariable(name, "f8", dims)
            v[:] = arr
            v.units = units
        ds.variables["somtp"].long_name = "slab-ocean temperature (docn.r somtp)"
        ds.source = str(Path(docn_r).resolve())
        ds.domain_file = str(Path(domain_file).resolve())
    return Path(out)


# ---------------------------------------------------------------------------
# cam_atm: cam.r temperature + water vapor (atmos.py)
# ---------------------------------------------------------------------------

PROFILE_SUFFIX = ".atmprofile.nc"
_CONST_RE = {
    "cpair": re.compile(r"^\s*CPDAIR:\s*([-+0-9.Ee]+)", re.M),
    "rair": re.compile(r"^\s*RAIR:\s*([-+0-9.Ee]+)", re.M),
    "zvir": re.compile(r"^\s*ZVIR:\s*([-+0-9.Ee]+)", re.M),
    "gravit": re.compile(r"SURFACE GRAVITY \(m/s2\):\s*([-+0-9.Ee]+)"),
    "ptop": re.compile(r"PTOP=\s*([-+0-9.Ee]+)"),
}


def locate_cam_restart(rundir) -> Path:
    """The cam.r file rpointer.atm names (its first line)."""
    rundir = Path(rundir)
    rp = rundir / "rpointer.atm"
    if not rp.is_file():
        raise FileNotFoundError(f"no rpointer.atm in {rundir}")
    lines = [l.strip() for l in rp.read_text().splitlines() if l.strip()]
    path = (rundir / lines[0]).resolve()
    if ".cam.r." not in path.name or not path.is_file():
        raise FileNotFoundError(f"rpointer.atm names {lines[0]}, not an existing cam.r")
    return path


def atm_constants_from_log(path):
    """AtmConstants from an ExoCAM atm.log (CPDAIR, RAIR, ZVIR, SURFACE
    GRAVITY, PTOP). ``path`` may be the log or a run directory (its newest
    atm.log.*)."""
    from .atmos import AtmConstants
    p = Path(path)
    if p.is_dir():
        logs = sorted(p.glob("atm.log.*"), key=lambda f: f.stat().st_mtime)
        logs = [f for f in logs if not f.name.endswith(".gz")]
        if not logs:
            raise FileNotFoundError(f"no atm.log.* in {p} (pass --atm-log)")
        p = logs[-1]
    text = p.read_text(errors="replace")
    vals = {}
    for k, rx in _CONST_RE.items():
        m = rx.search(text)
        if not m:
            raise ValueError(f"{p.name}: no {k} printout found (not an ExoCAM atm.log?)")
        vals[k] = float(m.group(1))
    return AtmConstants(**vals), p


def read_cam_state(path) -> Dict[str, np.ndarray]:
    from .atmos import STATE_FIELDS
    _require_netcdf()
    return {k: v.astype(float) for k, v in _read(path, list(STATE_FIELDS)).items()}


def apply_atmos_jump(path, dT_levels, constants, provenance: Optional[dict] = None,
                     dry_run: bool = False, profile_file=None,
                     extra_meta: Optional[dict] = None) -> JumpRecord:
    """Raise cam.r T by ``dT_levels`` (K per level) and q at fixed RH, in place
    with backup, keeping dry mass, TEOUT and the stratiform scheme's
    previous-step state consistent (``atmos.jump_state``). Same backup / log /
    re-apply-from-pristine rules as the other plugins; the record's factor slot
    holds the surface-level increment."""
    from .atmos import WRITTEN_FIELDS, CamAtmPlugin

    _require_netcdf()
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if ".cam.r." not in path.name:
        raise ValueError(f"{path.name} is not a cam.r file")
    backup, log = backup_path(path), log_path(path)
    if backup.exists():
        if is_jumped(backup):
            raise RuntimeError(f"{backup} is marked as already jumped; it is not "
                               f"a pristine backup — refusing")
        source = backup
    else:
        if is_jumped(path):
            raise RuntimeError(f"{path.name} was already jumped but its pristine "
                               f"backup {backup} is missing — refusing to shift twice")
        source = path

    before = read_cam_state(source)
    plugin = CamAtmPlugin(constants)
    jumped = plugin.apply_delta(before, np.asarray(dT_levels, dtype=float))
    after, report = plugin.enforce_constraints(before, jumped)
    rep = plugin.last_report

    date = restart_date(path)
    dT = np.asarray(dT_levels, dtype=float)
    meta = {
        "tool": "exocam-accelerate", "plugin": plugin.name,
        "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "case": case_of(path), "restart_file": path.name, "restart_date": date,
        "jump_model_year": first_model_year(date),
        "dT_surface_level": float(dT[-1]), "dT_levels": dT.tolist(),
        "fields": list(WRITTEN_FIELDS), "pristine_backup": str(backup),
        "constants": {k: getattr(constants, k) for k in
                      ("cpair", "rair", "zvir", "gravit", "ptop", "latvap", "latice")},
        "report": rep,
    }
    if extra_meta:
        meta.update(extra_meta)
    copy = None
    if profile_file is not None:
        copy = state_dir(path) / (path.name + PROFILE_SUFFIX)
        meta["profile"] = {"source": str(Path(profile_file).resolve()),
                           "sha256": file_sha256(profile_file), "copy": str(copy)}
    if dry_run:
        return JumpRecord(path, backup, log, float(dT[-1]), 1.0, False,
                          dict(report.adjustments), meta)

    backup.parent.mkdir(exist_ok=True)
    if copy is not None:
        shutil.copy2(profile_file, copy)
    if source is path:
        shutil.copy2(path, backup)
    else:
        shutil.copy2(backup, path)
    with netCDF4.Dataset(path, "r+") as ds:
        ds.set_auto_mask(False)
        for name in WRITTEN_FIELDS:
            var = ds.variables[name]
            var[:] = after[name].astype(var.dtype)
        ds.setncattr(ATTR, json.dumps({k: v for k, v in meta.items()
                                       if k not in ("dT_levels",)}))
    check = _read(path, list(WRITTEN_FIELDS))
    for name in WRITTEN_FIELDS:
        if not np.array_equal(check[name], after[name].astype(check[name].dtype)):
            raise RuntimeError(f"verification failed for {name} in {path}; restore "
                               f"with the pristine copy {backup}")
    log.write_text(json.dumps(dict(meta, adjustments=dict(report.adjustments),
                                   advice=provenance), indent=2))
    return JumpRecord(path, backup, log, float(dT[-1]), 1.0, True,
                      dict(report.adjustments), meta)


def write_profile_file(path, prof, dTS_jump: Optional[float] = None) -> str:
    _require_netcdf()
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lev", prof.gain.size)
        for name, arr, units in (("gain", prof.gain, "K/K"), ("raw_gain", prof.raw_gain, "K/K"),
                                 ("p_mid", prof.p_mid, "Pa"), ("T_now", prof.T_now, "K")):
            v = ds.createVariable(name, "f8", ("lev",))
            v[:] = arr
            v.units = units
        if prof.raw_gain_se is not None:
            v = ds.createVariable("raw_gain_se", "f8", ("lev",))
            v[:] = prof.raw_gain_se
            v.units = "K/K"
        if dTS_jump is not None:
            v = ds.createVariable("dT", "f8", ("lev",))
            v[:] = prof.gain * dTS_jump
            v.units = "K"
        ds.setncattr(ATTR, json.dumps(dict(prof.summary(), dTS_jump=dTS_jump)))
    return file_sha256(path)


def read_profile_file(path) -> Dict[str, np.ndarray]:
    _require_netcdf()
    return _read(path, ["gain", "raw_gain", "p_mid", "T_now"])
