"""Implementation of sampling algorithms."""
from __future__ import annotations

import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """A simple sampler that selects a random subset of structures."""

    def __init__(self, config: SamplingConfig) -> None:
        """
        Initialise the RandomSampler.

        Args:
            config: The Pydantic model for the sampling configuration.
        """
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Selects a random fraction of the input structures.

        Args:
            structures: The input list of `ase.Atoms` objects.

        Returns:
            A new list containing a randomly sampled subset of the input structures.
        """
        if self.config.method != "random":
            return structures

        num_to_select = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_select)
