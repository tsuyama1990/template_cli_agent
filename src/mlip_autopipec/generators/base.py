from abc import ABC, abstractmethod

from ase import Atoms

from mlip_autopipec.common.pydantic_models import SystemConfig


class BaseStructureGenerator(ABC):
    """Abstract base class for all structure generators."""

    def __init__(self, config: SystemConfig) -> None:
        self.config = config

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of initial atomic structures.

        This method must be implemented by all concrete generator classes.
        """
        raise NotImplementedError
