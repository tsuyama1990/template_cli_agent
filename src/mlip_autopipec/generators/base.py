"""Abstract base class for structure generators."""
from abc import ABC, abstractmethod
from typing import List

from ase.atoms import Atoms


class BaseStructureGenerator(ABC):
    """Abstract interface for all structure generator classes."""

    @abstractmethod
    def generate(self) -> List[Atoms]:
        """
        Generate a list of atomic structures.

        Returns:
            A list of ase.Atoms objects.
        """
        raise NotImplementedError
