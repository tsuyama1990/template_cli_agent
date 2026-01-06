"""Provides implementations of various sampling algorithms.

After the exploration stage generates a large number of candidate structures,
the sampling module is responsible for selecting a smaller, representative, and
information-rich subset. This module contains different samplers that can be
chosen via the user configuration. For Cycle 1, only a simple `RandomSampler`
is implemented.
"""

from __future__ import annotations

import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """Selects a random subset of the provided atomic structures.

    This sampler implements a basic random sampling strategy. It takes a list of
    structures and returns a new list containing a randomly chosen fraction of
    the original structures, as specified by the user's configuration.

    Attributes:
        config (SamplingConfig): The validated Pydantic model for the
            sampling stage configuration.

    """

    def __init__(self, config: SamplingConfig) -> None:
        """Initialise the RandomSampler.

        Args:
            config: The Pydantic model for the sampling configuration.

        """
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """Select a random fraction of the input structures.

        Args:
            structures: The input list of `ase.Atoms` objects.

        Returns:
            A new list containing a randomly sampled subset of the input
            structures. Returns the original list if the method is not 'random'.

        """
        if self.config.method != "random":
            return structures

        num_to_select = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_select)
