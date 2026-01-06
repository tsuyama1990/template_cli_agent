"""
Defines the base class for all structure generators.
"""

from abc import ABC, abstractmethod

from ase import Atoms


class BaseStructureGenerator(ABC):
    """
    Abstract base class for all structure generators.

    Each generator must implement the `generate` method, which returns a list
    of initial atomic structures (ASE Atoms objects).
    """

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of initial atomic structures.

        Returns:
            A list of ASE Atoms objects representing the generated structures.
        """
