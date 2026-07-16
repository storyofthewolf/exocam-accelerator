"""exocam-accelerate: safeguarded convergence acceleration for ExoCAM/CESM.

Core operation: X_new = X + <dX/dt> * Δt (Wordsworth et al. 2013 §2.3;
Turbet et al. 2021 Methods §4), wrapped in safeguards that decide when the
extrapolation is trustworthy and how large it may be. See README.md.
"""

from .hindcast import (
    GateCalibration,
    HindcastResult,
    annual_mean_series,
    calibrate_gate,
    default_origins,
    error_by_dt,
    hindcast_at,
    max_safe_dt,
    sweep_hindcasts,
)
from .safeguards import ClipReport, GateConfig, GateResult, assess_trustworthiness, clip_step
from .schedule import DtSchedule
from .stepper import StepProposal, propose_step
from .trend_io import load_case, read_trend_text, trend_series
from .trends import (
    TrendSeries,
    build_area_weights,
    fit_curvature,
    fit_tendency,
    horizontal_mean,
    layer_mean_series,
    segment_slopes,
)

__version__ = "0.1.0"

__all__ = [
    "ClipReport",
    "DtSchedule",
    "GateCalibration",
    "GateConfig",
    "GateResult",
    "HindcastResult",
    "StepProposal",
    "TrendSeries",
    "annual_mean_series",
    "assess_trustworthiness",
    "build_area_weights",
    "calibrate_gate",
    "clip_step",
    "default_origins",
    "error_by_dt",
    "fit_curvature",
    "fit_tendency",
    "hindcast_at",
    "horizontal_mean",
    "layer_mean_series",
    "load_case",
    "max_safe_dt",
    "propose_step",
    "read_trend_text",
    "segment_slopes",
    "sweep_hindcasts",
    "trend_series",
]
