"""Variable-plugin interface for post-jump physical-constraint re-imposition.

Implemented plugins: ``aqua_ice`` (``aqua_ice.py``), for the in-place
continuation path decided 2026-09-25 (docs/restart-integration-questions.md
§7). Other targets (somtp, clm) remain design-only.

Safeguard 5 of the design: after an extrapolation jump, each accelerated
variable must have its physical constraints re-imposed. This is the only
variable-specific piece of the whole scheme. Precedent (Wordsworth et al.
2013): after each ice jump, remove ice where the extrapolated thickness went
negative or coverage was seasonal, then renormalize so total water mass is
conserved.

A plugin owns everything variable-specific:
  * which restart file type its fields live in (cam.r, cice.r, clm2.r, ...),
  * which variables it reads and perturbs,
  * how a per-layer (or pointwise) delta maps onto the gridded fields,
  * which physical constraints to re-impose after the jump.

The first planned plugin is aqua-ice (snowball sea-ice energy/volume:
``eicen``, ``vicen``, ``esnon``, ``vsnon`` in cice.r), replacing the
prototype's bare multiplicative scaling.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple

import numpy as np


@dataclass(frozen=True)
class ConstraintReport:
    """What a plugin's constraint pass changed, for the step log."""

    plugin: str
    adjustments: Mapping[str, str]  # variable name -> human-readable summary


class VariablePlugin(ABC):
    """One acceleration target: its fields, its delta mapping, its physics.

    Implementations operate on in-memory arrays only. File reading/writing
    stays in the (future) restart layer, which will call these hooks.
    """

    #: short identifier used in configs and logs, e.g. "aqua_ice"
    name: str

    #: restart file type the fields live in, e.g. "cice.r"
    file_kind: str

    #: names of the restart variables this plugin perturbs
    variables: Tuple[str, ...]

    @abstractmethod
    def apply_delta(
        self,
        fields: Dict[str, np.ndarray],
        delta: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Map an accepted StepProposal delta onto the gridded fields.

        ``fields`` maps each name in ``variables`` to its array as read from
        the restart file; ``delta`` is the stepper's clipped increment (e.g.
        per-layer, to be broadcast across the horizontal axes; or pointwise
        for Wordsworth-style surface fields). Returns new arrays; must not
        mutate the inputs.
        """
        raise NotImplementedError

    @abstractmethod
    def enforce_constraints(
        self,
        fields_before: Dict[str, np.ndarray],
        fields_after: Dict[str, np.ndarray],
    ) -> Tuple[Dict[str, np.ndarray], ConstraintReport]:
        """Re-impose physical constraints on post-jump fields.

        Receives the pre-jump fields (for conserved-quantity budgets) and the
        jumped fields; returns corrected fields plus a report of what was
        adjusted. Example (aqua-ice): clamp ice energy/volume at physical
        bounds, keep energy-volume consistency, conserve total water mass.
        """
        raise NotImplementedError


#: Registry mapping plugin names to instances. Populated once real plugins
#: exist; the acceleration driver will look targets up here.
PLUGIN_REGISTRY: Dict[str, VariablePlugin] = {}
