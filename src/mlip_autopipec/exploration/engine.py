"""Placeholder for the MD Exploration Engine."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import ExplorationConfig


class MDExplorer:
    """A placeholder for the Molecular Dynamics exploration engine."""

    def __init__(self, config: ExplorationConfig) -> None:
        """
        Initialise the MDExplorer.

        Args:
            config: The Pydantic model for the exploration configuration.
        """
        self.config = config

    def run_md(self, structures: list[Atoms]) -> list[Atoms]:
        """
        A placeholder method that will eventually run MD simulations.

        In Cycle 1, this method is a simple pass-through.

        Args:
            structures: The input list of `ase.Atoms` objects.

        Returns:
            The same list of `ase.Atoms` objects, unmodified.
        """
        print("Exploration stage: Placeholder - passing structures through.")
        return structures
