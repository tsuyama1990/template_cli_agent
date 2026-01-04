from abc import ABC, abstractmethod

from ase import Atoms


class BaseStructureGenerator(ABC):
    """
    Abstract base class for all structure generators.
    """

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of atomic structures.
        """
