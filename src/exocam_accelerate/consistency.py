"""Post-step consistency-check interface — safeguard 6.

DESIGN ONLY — the run-forward half does not exist yet, so this module defines
the contract and nothing else.

The idea: an accepted extrapolation step is a *prediction* that the model,
restarted from the jumped state, will keep drifting along the same trajectory
(same tendency sign, comparable magnitude, from the new value). If instead the
model immediately corrects *away* from the jump — tendency reversed, or the
state relaxing back toward the pre-jump value — the step crossed a
nonlinearity that the pre-step gate could not see, and the jump must be
rolled back to the pre-jump restart set.

Intended flow (future driver):

    proposal = propose_step(...)                    # implemented (stepper)
    record   = PredictionRecord.from_proposal(...)  # implemented below
    ... apply jump, run model forward, measure new trend window ...
    verdict  = checker.evaluate(record, post_series)  # NOT implemented
    if verdict.action is Action.ROLLBACK: restore backup restart set
"""

from __future__ import annotations

import enum
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .stepper import StepProposal
from .trends import TrendSeries


@dataclass(frozen=True)
class PredictionRecord:
    """What an accepted step predicted, kept for post-run comparison."""

    variable: str
    dt_years: float
    pre_step_value: np.ndarray   # state (e.g. per-layer means) before the jump
    tendency: np.ndarray         # fitted <dX/dt> the jump was built from
    delta: np.ndarray            # the clipped increment actually applied

    @classmethod
    def from_proposal(cls, variable: str, series: TrendSeries,
                      proposal: StepProposal) -> "PredictionRecord":
        if not proposal.accepted:
            raise ValueError("cannot record a prediction for a refused step")
        return cls(
            variable=variable,
            dt_years=proposal.dt_years,
            pre_step_value=np.asarray(series.values[-1]),
            tendency=proposal.tendency,
            delta=proposal.delta,
        )


class Action(enum.Enum):
    ACCEPT = "accept"      # model continued along the predicted trajectory
    ROLLBACK = "rollback"  # model corrected away from the jump: restore backup


@dataclass(frozen=True)
class ConsistencyVerdict:
    action: Action
    #: fraction-of-prediction agreement metric, definition TBD with the
    #: implementation (open question: threshold and per-layer vs aggregate)
    score: Optional[float] = None
    reasons: tuple = ()


class ConsistencyCheck(ABC):
    """Contract for the post-step consistency evaluation.

    ``post_series`` is a trend window measured from the forward run *after*
    the jump, built the same way as the pre-step window (per-layer horizontal
    means for the default pipeline). Implementations decide how to compare it
    against the prediction — at minimum the post-run tendency must not have
    reversed sign against ``record.tendency`` while the state moves back
    toward ``record.pre_step_value``.
    """

    @abstractmethod
    def evaluate(self, record: PredictionRecord,
                 post_series: TrendSeries) -> ConsistencyVerdict:
        raise NotImplementedError
