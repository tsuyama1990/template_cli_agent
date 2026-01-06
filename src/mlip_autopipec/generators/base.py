from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ase import Atoms


class BaseStructureGenerator(ABC):
    """Abstract base class for all structure generators."""

    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.elements: list[str] = config.get("elements", [])
        self.num_structures: int = config.get("num_structures", 1)

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of atomic structures.

        This method must be implemented by all concrete generator classes.

        Returns:
            A list of `ase.Atoms` objects.
        """
        raise NotImplementedError
