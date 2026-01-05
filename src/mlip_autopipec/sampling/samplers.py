"""Implementation of sampling algorithms."""

import random

from ase.atoms import Atoms

from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """A sampler that selects a random subset of structures."""

    def __init__(self, config: SamplingConfig) -> None:
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Selects a random subset of the input structures.

        Args:
            structures: A list of ASE Atoms objects.

        Returns:
            A new list of ASE Atoms objects containing a random subset of the input.
        """
        num_to_sample = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_sample)
