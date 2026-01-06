"""Abstract base class for structure generators."""
from __future__ import annotations

import abc
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import SystemConfig


class BaseStructureGenerator(abc.ABC):
    """Abstract base class for all structure generators."""

    def __init__(self, config: SystemConfig) -> None:
        """
        Initialise the generator with the system configuration.

        Args:
            config: The Pydantic model containing the system configuration.
        """
        self.config = config

    @abc.abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of atomic structures.

        This method must be implemented by all concrete generator classes.

        Returns:
            A list of `ase.Atoms` objects.
        """
        raise NotImplementedError
