# -*- coding: utf-8 -*-
"""
Placeholder MDExplorer for the exploration stage.

This module provides a simple, pass-through implementation of the MDExplorer
for Cycle 1, fulfilling the architectural requirements without implementing
the full complexity of MD simulations.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

logging.basicConfig(level=logging.INFO)


class MDExplorer:
    """
    Placeholder for the MD exploration engine.

    In Cycle 1, this class serves as a simple pass-through component in the
    pipeline. Its `run_md` method returns the input structures without
    modification, allowing the end-to-end pipeline to be tested.
    """

    def run_md(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Run the (placeholder) MD simulation.

        This method logs a message and returns the input list of structures
        unmodified, acting as a placeholder for the real MD simulation.

        Args:
            structures: The list of input structures.

        Returns:
            The list of structures, unmodified.
        """
        logging.info("Running placeholder MDExplorer: returning structures as is.")
        return structures
