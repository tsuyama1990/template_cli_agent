"""Placeholder for the Molecular Dynamics (MD) Exploration Engine.

This module contains the `MDExplorer` class, which will eventually be responsible
for running MD simulations to explore the potential energy surface of the
material. In Cycle 1, this component acts as a simple pass-through placeholder
to allow for the end-to-end testing of the pipeline's architecture. Its
interface is defined with the final design in mind for forward compatibility.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import ExplorationConfig


class MDExplorer:
    """A placeholder for the Molecular Dynamics exploration engine.

    In Cycle 1, this class does not perform any simulations. It simply accepts
    a list of structures and returns it unmodified, logging a message to
    indicate that it was called.

    Attributes:
        config (ExplorationConfig): The validated Pydantic model for the
            exploration stage configuration.

    """

    def __init__(self, config: ExplorationConfig) -> None:
        """Initialise the MDExplorer.

        Args:
            config: The Pydantic model for the exploration configuration.

        """
        self.config = config

    def run_md(self, structures: list[Atoms]) -> list[Atoms]:
        """Run the (placeholder) MD simulation.

        This method will, in a future cycle, execute MD runs for each of the
        input structures. For now, it is a simple pass-through.

        Args:
            structures: The input list of `ase.Atoms` objects.

        Returns:
            The same list of `ase.Atoms` objects, unmodified.

        """
        logging.info("Exploration stage: Placeholder - passing structures through.")
        return structures
