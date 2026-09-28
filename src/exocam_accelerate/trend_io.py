"""Adapter for exocam-trend time-series output — the only file I/O in the core.

exocam-trend (sibling repo, a dependency) writes global-mean time series as
whitespace-separated text with a header row of ``month`` followed by
``VAR_native  VAR_int1  VAR_int2`` triplets per variable (native monthly
value, short-window running mean, long-window running mean). This module
parses that format into TrendSeries objects; it does not read model netCDF
output itself.

exocam-trend currently produces only global means of 2D fields — no per-layer
(vertical) differentiation. Per-layer trend series for 3D acceleration targets
must be built directly from model output via trends.layer_mean_series until
exocam-trend grows that capability.
"""

from __future__ import annotations

import glob
import hashlib
import re
from pathlib import Path
from typing import Dict, List

import numpy as np

from .trends import TrendSeries

MONTHS_PER_YEAR = 12.0

# exocam-trend writes one file per component per case, suffixed by stream.
# A case is the set of these files sharing a <case>_<first>-<last> stem.
_COMPONENT_SUFFIXES = ("cam", "cice", "clm")


def read_trend_text(path) -> Dict[str, np.ndarray]:
    """Parse an exocam-trend data/*.txt file into named columns.

    Returns a dict mapping each header name (``month``, ``TS_native``,
    ``TS_int1``, ...) to a 1D float array.
    """
    path = Path(path)
    with path.open() as f:
        header: List[str] = f.readline().split()
        if not header or header[0] != "month":
            raise ValueError(
                f"{path} does not look like exocam-trend output "
                f"(first header column is {header[:1]!r}, expected 'month')"
            )
        data = np.loadtxt(f, ndmin=2)
    if data.shape[1] != len(header):
        raise ValueError(
            f"{path}: {data.shape[1]} data columns != {len(header)} header names"
        )
    return {name: data[:, i] for i, name in enumerate(header)}


def trend_series(columns: Dict[str, np.ndarray], variable: str,
                 which: str = "native", last_n_months: int = None) -> TrendSeries:
    """Build a TrendSeries for one variable from parsed exocam-trend columns.

    which selects the column flavor: 'native' (monthly instantaneous),
    'int1' (short running mean), or 'int2' (long running mean). Running-mean
    columns lag the true state early in the series; prefer 'native' over a
    deliberately chosen trend window (``last_n_months``), or slice upstream.

    Times are converted from month index to years.
    """
    key = f"{variable}_{which}"
    if key not in columns:
        available = sorted(k for k in columns if k != "month")
        raise KeyError(f"column {key!r} not found; available: {available}")
    months = columns["month"]
    values = columns[key]
    if last_n_months is not None:
        months = months[-last_n_months:]
        values = values[-last_n_months:]
    return TrendSeries(times=months / MONTHS_PER_YEAR, values=values)


_SPAN = re.compile(r"_(\d{4})-(\d{2})-(\d{4})-(\d{2})_(?:cam|cice|clm)\.txt$")


def case_start_year(directory, case_id: str) -> int:
    """Model year of the first month in a case's exocam-trend files.

    Parsed from the ``<case>_<YYYY>-<MM>-<YYYY>-<MM>_<comp>.txt`` file name.
    Annual means are binned from the series' first month, so they align with
    model years only when the series starts in January; anything else raises.
    """
    for path in sorted(Path(directory).glob(f"{case_id}_*.txt")):
        m = _SPAN.search(path.name)
        if m:
            if m.group(2) != "01":
                raise ValueError(f"{path.name}: series starts in month {m.group(2)}; "
                                 f"annual means need a January start")
            return int(m.group(1))
    raise FileNotFoundError(f"no exocam-trend files for {case_id} in {directory}")


def file_provenance(directory, case_id: str) -> Dict[str, str]:
    """``{file name: sha256 hex digest}`` for a case's exocam-trend files.

    Binds advice (``advise.Advice.to_dict()["provenance"]``) to the exact
    input series it was fit from, for later audit (feasibility-review, Stage
    0 item 5). Uses the same file set as ``load_case``; returns ``{}`` if
    none are found (``load_case`` itself will raise on that).
    """
    directory = Path(directory)
    out: Dict[str, str] = {}
    for suffix in _COMPONENT_SUFFIXES:
        for path in sorted(directory.glob(f"{case_id}_*_{suffix}.txt")):
            out[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def load_case(directory, case_id: str) -> Dict[str, np.ndarray]:
    """Merge a case's per-component exocam-trend files into one column dict.

    exocam-trend writes ``<case>_<first>-<last>_{cam,cice,clm}.txt`` (the
    sea-ice/snowpack variables live in the ``_cice`` file, ICEFRAC/energy in
    ``_cam``, land in ``_clm``). This finds every component file whose name
    starts with ``case_id`` in ``directory`` and returns their columns merged
    under a shared ``month`` axis, so a variable can be requested without the
    caller knowing which component emitted it.

    Raises if the component files disagree on their month axis (they should be
    identical when generated together by ``run_trend_batch.sh``).
    """
    directory = Path(directory)
    paths = [p for suffix in _COMPONENT_SUFFIXES
             for p in sorted(glob.glob(str(directory / f"{case_id}_*_{suffix}.txt")))]
    if not paths:
        raise FileNotFoundError(
            f"no exocam-trend files matching {case_id}_*_{{cam,cice,clm}}.txt "
            f"in {directory}"
        )
    return merge_trend_files(paths, case_id)


def merge_trend_files(paths, case_id: str = "?") -> Dict[str, np.ndarray]:
    """Merge explicit per-component exocam-trend files under one ``month`` axis.

    The worker behind ``load_case``; also used directly when a directory holds
    several spans of one case and the caller has picked which to read.
    """
    merged: Dict[str, np.ndarray] = {}
    month_axis = None
    for path in paths:
        cols = read_trend_text(path)
        if month_axis is None:
            month_axis = cols["month"]
            merged["month"] = month_axis
        elif not np.array_equal(cols["month"], month_axis):
            raise ValueError(
                f"{path}: month axis differs from earlier component file(s) "
                f"for case {case_id!r}; were they generated in one run?"
            )
        for name, arr in cols.items():
            if name == "month":
                continue
            if name in merged:
                raise ValueError(
                    f"{path}: column {name!r} already loaded for case "
                    f"{case_id!r} from another component"
                )
            merged[name] = arr
    if month_axis is None:
        raise FileNotFoundError(f"no exocam-trend files given for case {case_id!r}")
    return merged
