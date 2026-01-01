"""Defines the base class for all structure generators."""
from abc import ABC, abstractmethod

from ase import Atoms

from mlip_autopipec.config import FullConfig


class BaseStructureGenerator(ABC):
    """
    Abstract base class for all structure generation modules.

    Each generator is responsible for creating a list of atomic structures
    based on the provided configuration. The logic is self-contained, and
    the class does not interact with any external systems like databases.
    """

    def __init__(self, config: FullConfig):
        """
        Initializes the base generator.

        Args:
            config: The full, expanded configuration object.
        """
        self.config = config
        self.system_info = config.system

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        The main method to generate a list of atomic structures.

        This method must be implemented by all subclasses.

        Returns:
            A list of ase.Atoms objects.
        """
        pass
