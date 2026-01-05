"""Abstract base class for structure generators."""

from abc import ABC, abstractmethod

from ase.atoms import Atoms


class BaseStructureGenerator(ABC):
    """Abstract base class for structure generators."""

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """Generate a list of atomic structures."""
