"""Defines the abstract interface for all structure generators.

This module provides the `BaseStructureGenerator`, an abstract base class (ABC)
that establishes a common contract for all components responsible for creating
the initial "seed" atomic structures. This polymorphic design is a key part of
the application's extensible architecture, allowing new types of structure
generators (e.g., for ionic crystals, surfaces, or molecules) to be added
without modifying the core pipeline orchestration logic.
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import SystemConfig


class BaseStructureGenerator(abc.ABC):
    """Abstract base class for all structure generators.

    This class defines the interface that all concrete generator classes must
    implement. It ensures that the main `PipelineRunner` can work with any
    type of generator through a consistent API.

    Attributes:
        config (SystemConfig): The validated Pydantic model containing the
            configuration for the physical system.

    """

    def __init__(self, config: SystemConfig) -> None:
        """Initialise the generator with the system configuration.

        Args:
            config: The Pydantic model containing the system configuration.

        """
        self.config = config

    @abc.abstractmethod
    def generate(self) -> list[Atoms]:
        """Generate a list of atomic structures.

        This method must be implemented by all concrete generator classes.

        Returns:
            A list of `ase.Atoms` objects.

        """
        raise NotImplementedError
