"""Decreasing Δt schedules for the acceleration loop.

Both precedents ramp the acceleration timestep down as convergence
approaches — large jumps early when far from equilibrium, small jumps near
it. Wordsworth et al. 2013 used Δt = 100 yr for the first 5 iterations, then
10 yr for the final 15; Turbet et al. 2021 stepped n_days through 50 → 10 → 0.
A trailing Δt of 0 is allowed and means "stop accelerating, run normally".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple


@dataclass(frozen=True)
class DtSchedule:
    """A non-increasing sequence of (n_iterations, dt_years) phases."""

    phases: Tuple[Tuple[int, float], ...]

    def __post_init__(self):
        phases = tuple((int(n), float(dt)) for n, dt in self.phases)
        if not phases:
            raise ValueError("schedule needs at least one phase")
        for n, dt in phases:
            if n < 1:
                raise ValueError(f"phase iteration count must be >= 1, got {n}")
            if dt < 0:
                raise ValueError(f"dt must be >= 0, got {dt}")
        dts = [dt for _, dt in phases]
        if any(later > earlier for earlier, later in zip(dts, dts[1:])):
            raise ValueError(
                f"dt must be non-increasing across phases, got {dts}: "
                "the whole point is smaller jumps near equilibrium"
            )
        if any(dt == 0 for _, dt in phases[:-1]):
            raise ValueError("only the final phase may have dt = 0")
        object.__setattr__(self, "phases", phases)

    @classmethod
    def from_phases(cls, phases: Sequence[Tuple[int, float]]) -> "DtSchedule":
        return cls(tuple(phases))

    @classmethod
    def wordsworth2013(cls) -> "DtSchedule":
        """The published early-Mars ice schedule: 5 × 100 yr, then 15 × 10 yr."""
        return cls(((5, 100.0), (15, 10.0)))

    @property
    def total_iterations(self) -> int:
        return sum(n for n, _ in self.phases)

    def dt_for(self, iteration: int) -> float:
        """Δt (years) for a 0-based iteration index.

        Raises IndexError past the end of the schedule — the loop is over;
        callers must not invent further steps.
        """
        if iteration < 0:
            raise IndexError("iteration must be >= 0")
        remaining = iteration
        for n, dt in self.phases:
            if remaining < n:
                return dt
            remaining -= n
        raise IndexError(
            f"iteration {iteration} beyond schedule end ({self.total_iterations})"
        )

    def __iter__(self):
        for n, dt in self.phases:
            for _ in range(n):
                yield dt
