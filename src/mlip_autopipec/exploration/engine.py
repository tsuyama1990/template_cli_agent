from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import ExplorationConfig

logger = logging.getLogger(__name__)


class MDExplorer:
    """A placeholder for the molecular dynamics exploration engine."""

    def __init__(self, config: ExplorationConfig):
        """
        Initialise the MDExplorer.

        Args:
            config: The exploration configuration.
        """
        self.config = config

    def explore(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Run the (placeholder) exploration.

        Args:
            structures: The input structures.

        Returns:
            The structures, unmodified.
        """
        logger.info("Skipping MD exploration in Cycle 1.")
        return structures
