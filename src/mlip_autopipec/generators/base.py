# src/mlip_autopipec/generators/base.py
"""
Defines the abstract base class for all structure generators, enforcing a common
interface and providing shared validation logic.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig

if TYPE_CHECKING:
    pass


class BaseStructureGenerator(ABC):
    """Abstract base class for structure generators."""

    def __init__(self, config: FullConfig) -> None:
        self.config = config

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """Generate a list of atomic structures."""
        raise NotImplementedError

    def _validate_structure(self, atoms: Atoms) -> bool:
        """
        Perform physics-based validation checks on a generated structure.
        For Cycle 1, this includes a simple overlap check.
        """
        # Simple overlap check: no two atoms should be closer than 1.0 Angstrom
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if atoms.get_distance(i, j) < 1.0:  # type: ignore[no-untyped-call]
                    return False
        return True
