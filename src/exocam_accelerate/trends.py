"""Trend measurement: area weights, per-layer horizontal means, tendency fits.

Tendencies feeding an acceleration step are measured from horizontally
averaged per-layer trends (Turbet et al. 2021, Methods §4), not pointwise
gridcell trends: pointwise tendencies amplify dynamical noise, while the
horizontal mean isolates the slow coherent drift toward equilibrium.

The fitting helpers are shape-generic over trailing axes — values of shape
(nt,), (nt, nlev), or (nt, nlat, nlon) all work — so variable plugins that
genuinely need pointwise trends (Wordsworth's surface-ice case) can reuse
them. The default pipeline is per-layer means.

Time is in years everywhere; fitted tendencies are per year.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def build_area_weights(lon, lat):
    """Normalized 2D area weights for a regular lat/lon grid.

    Adapted from exocam-trend's ``trend_core.build_area_weights`` (Wolf),
    duplicated here because exocam-trend is not an importable package.

    Returns an array of shape (nlat, nlon) summing to 1.0. Cell edges sit at
    latitude midpoints, with poles at ±90°.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lat_rad = np.deg2rad(lat)

    lat_edges = np.empty(len(lat) + 1)
    lat_edges[0] = -np.pi / 2.0
    lat_edges[-1] = np.pi / 2.0
    lat_edges[1:-1] = 0.5 * (lat_rad[:-1] + lat_rad[1:])

    dlon_rad = np.deg2rad(lon[1] - lon[0]) if len(lon) > 1 else 2.0 * np.pi

    d_sin_lat = np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1])
    weights = d_sin_lat[:, np.newaxis] * np.full(len(lon), dlon_rad)

    return weights / weights.sum()


def horizontal_mean(field, weights):
    """Area-weighted horizontal mean over the trailing (nlat, nlon) axes.

    ``field`` may carry any leading axes — (nlat, nlon), (nlev, nlat, nlon),
    (nt, nlev, nlat, nlon) — and the result drops the trailing two. Masked
    arrays are honored: masked cells are excluded from numerator and
    denominator alike.
    """
    field = np.ma.asarray(field)
    if field.shape[-2:] != weights.shape:
        raise ValueError(
            f"trailing axes {field.shape[-2:]} do not match weights {weights.shape}"
        )
    w = np.broadcast_to(weights, field.shape)
    masked = np.ma.getmaskarray(field)
    wsum = np.where(masked, 0.0, w).sum(axis=(-2, -1))
    total = (field.filled(0.0) * np.where(masked, 0.0, w)).sum(axis=(-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / wsum
    return np.where(wsum > 0, out, np.nan)


@dataclass(frozen=True)
class TrendSeries:
    """A trend window: sample times (years) and values at each time.

    ``values`` has shape (nt, *state_shape) — state_shape is () for a global
    scalar, (nlev,) for per-layer means, or a grid shape for pointwise use.
    """

    times: np.ndarray
    values: np.ndarray

    def __post_init__(self):
        times = np.asarray(self.times, dtype=float)
        values = np.asarray(self.values, dtype=float)
        if times.ndim != 1:
            raise ValueError("times must be 1D")
        if len(times) < 2:
            raise ValueError("a trend window needs at least 2 samples")
        if values.shape[0] != len(times):
            raise ValueError(
                f"values first axis {values.shape[0]} != len(times) {len(times)}"
            )
        if np.any(np.diff(times) <= 0):
            raise ValueError("times must be strictly increasing")
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "values", values)

    @property
    def window_length(self) -> float:
        return float(self.times[-1] - self.times[0])

    @property
    def state_shape(self):
        return self.values.shape[1:]


def layer_mean_series(times, fields, weights) -> TrendSeries:
    """Build a per-layer-mean TrendSeries from gridded snapshots.

    ``fields`` has shape (nt, nlev, nlat, nlon) or (nt, nlat, nlon); the
    horizontal axes are averaged away, leaving (nt, nlev) or (nt,).
    """
    return TrendSeries(times=np.asarray(times, dtype=float),
                       values=horizontal_mean(fields, weights))


def _polyfit(times, values, degree):
    """Least-squares polynomial fit along axis 0, generic over trailing axes.

    Times are centered on their midpoint before fitting for numerical
    conditioning. Returns coefficients (degree+1, *state_shape), lowest order
    first, valid in the centered time coordinate.
    """
    t = times - 0.5 * (times[0] + times[-1])
    vander = np.vander(t, degree + 1, increasing=True)
    flat = values.reshape(len(times), -1)
    coeffs, *_ = np.linalg.lstsq(vander, flat, rcond=None)
    return coeffs.reshape(degree + 1, *values.shape[1:])


def fit_tendency(series: TrendSeries):
    """Linear least-squares tendency <dX/dt> per trailing element (per year)."""
    coeffs = _polyfit(series.times, series.values, 1)
    return coeffs[1]


def fit_curvature(series: TrendSeries):
    """Quadratic fit; returns (tendency b, curvature d²X/dt² = 2c).

    b is the slope at the window midpoint; 2c is the constant second
    derivative of the fitted parabola. Needs at least 3 samples.
    """
    if len(series.times) < 3:
        raise ValueError("curvature fit needs at least 3 samples")
    coeffs = _polyfit(series.times, series.values, 2)
    return coeffs[1], 2.0 * coeffs[2]


def segment_slopes(series: TrendSeries, n_segments: int):
    """Linear slopes over ``n_segments`` contiguous sub-windows.

    Returns shape (n_segments, *state_shape). Used by the trustworthiness
    gate's sign-stability check. Every segment must contain ≥ 2 samples.
    """
    nt = len(series.times)
    if n_segments < 2:
        raise ValueError("need at least 2 segments")
    if nt < 2 * n_segments:
        raise ValueError(
            f"{nt} samples cannot form {n_segments} segments of >= 2 samples"
        )
    bounds = np.linspace(0, nt, n_segments + 1).astype(int)
    slopes = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        seg = _polyfit(series.times[lo:hi], series.values[lo:hi], 1)
        slopes.append(seg[1])
    return np.stack(slopes)
