from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


class BaseStructureGenerator(ABC):
    """Abstract base class for structure generators."""

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """Generates a list of ASE Atoms objects."""
        raise NotImplementedError
