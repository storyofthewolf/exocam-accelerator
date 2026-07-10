"""exocam-accelerate: safeguarded convergence acceleration for ExoCAM/CESM.

Core operation: X_new = X + <dX/dt> * Δt (Wordsworth et al. 2013 §2.3;
Turbet et al. 2021 Methods §4), wrapped in safeguards that decide when the
extrapolation is trustworthy and how large it may be. See README.md.
"""

from .safeguards import ClipReport, GateConfig, GateResult, assess_trustworthiness, clip_step
from .schedule import DtSchedule
from .stepper import StepProposal, propose_step
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
    "GateConfig",
    "GateResult",
    "StepProposal",
    "TrendSeries",
    "assess_trustworthiness",
    "build_area_weights",
    "clip_step",
    "fit_curvature",
    "fit_tendency",
    "horizontal_mean",
    "layer_mean_series",
    "propose_step",
    "segment_slopes",
]
