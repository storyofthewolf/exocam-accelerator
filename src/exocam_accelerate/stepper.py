"""The safeguarded extrapolation pipeline: gate → fit → extrapolate → clip.

Pure array computation. The result is a proposed per-element increment
(per-layer for the default horizontal-mean pipeline); applying it to gridded
restart fields, re-imposing physical constraints (plugins), and writing files
belong to later layers that do not exist yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .safeguards import ClipReport, GateConfig, GateResult, assess_trustworthiness, clip_step
from .trends import TrendSeries, fit_tendency


@dataclass(frozen=True)
class StepProposal:
    """Outcome of a proposed acceleration step for one variable.

    accepted is False when the trustworthiness gate refused; then delta is
    None and gate.reasons says why. When accepted, ``delta`` has the series'
    state shape (e.g. per-layer) and is already hard-clipped; ``clip`` reports
    whether and where the limit bit.
    """

    accepted: bool
    gate: GateResult
    dt_years: float
    tendency: Optional[np.ndarray] = None
    delta: Optional[np.ndarray] = None
    clip: Optional[ClipReport] = None

    def __bool__(self) -> bool:
        return self.accepted


def propose_step(
    series: TrendSeries,
    dt_years: float,
    max_abs_step: float,
    gate_config: GateConfig = GateConfig(),
) -> StepProposal:
    """Propose X_new − X = clip(<dX/dt> · Δt) for one variable, or refuse.

    ``series`` is normally the per-layer horizontal-mean trend window
    (trends.layer_mean_series). A refusal means the trend window sits too
    close to a nonlinear feedback threshold for linear extrapolation; the
    correct response upstream is to keep running the model normally, not to
    retry with a smaller Δt.
    """
    if dt_years < 0:
        raise ValueError("dt_years must be >= 0")

    gate = assess_trustworthiness(series, gate_config)
    if not gate:
        return StepProposal(accepted=False, gate=gate, dt_years=float(dt_years))

    tendency = fit_tendency(series)
    report = clip_step(tendency * dt_years, max_abs_step)
    return StepProposal(
        accepted=True,
        gate=gate,
        dt_years=float(dt_years),
        tendency=np.asarray(tendency),
        delta=report.clipped,
        clip=report,
    )
