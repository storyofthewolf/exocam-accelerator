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

from pathlib import Path
from typing import Dict, List

import numpy as np

from .trends import TrendSeries

MONTHS_PER_YEAR = 12.0


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
