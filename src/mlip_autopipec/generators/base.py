"""Abstract base class for structure generators."""
from abc import ABC, abstractmethod
from typing import List

from ase.atoms import Atoms


class BaseStructureGenerator(ABC):
    """Defines the interface for all structure generators."""

    @abstractmethod
    def generate(self) -> List[Atoms]:
        """
        Generate a list of atomic structures.

        Returns:
            A list of ase.Atoms objects.
        """
        raise NotImplementedError
