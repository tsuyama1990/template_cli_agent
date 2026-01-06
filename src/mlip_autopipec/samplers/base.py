# src/mlip_autopipec/samplers/base.py
"""Abstract base class for all sampling algorithms."""

from abc import ABC, abstractmethod

from ase import Atoms


class BaseSampler(ABC):
    """Defines the interface for all sampling algorithm implementations."""

    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Selects a subset of structures from a given trajectory.

        Args:
            trajectory: A list of ASE Atoms objects from the exploration phase.

        Returns:
            A list of selected ASE Atoms objects.
        """
