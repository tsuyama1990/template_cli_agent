# src/mlip_autopipec/generators/base.py
"""Abstract base class for all structure generators."""

from abc import ABC, abstractmethod

from ase import Atoms


class BaseStructureGenerator(ABC):
    """Defines the interface for all structure generator implementations."""

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generates a list of initial atomic structures.

        Returns:
            A list of ASE Atoms objects representing the generated structures.
        """
