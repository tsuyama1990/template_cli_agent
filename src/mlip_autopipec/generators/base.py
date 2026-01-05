"""Abstract BaseStructureGenerator class, defining the interface."""
from abc import ABC, abstractmethod

from ase import Atoms


class BaseStructureGenerator(ABC):
    """Abstract base class for structure generators."""

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """Generate a list of structures."""
        pass
