"""Implementation of sampling algorithms."""
import random
from typing import List
from ase.atoms import Atoms
from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """Selects a random subset of structures."""

    def __init__(self, config: SamplingConfig):
        self.config = config

    def sample(self, structures: List[Atoms]) -> List[Atoms]:
        """
        Randomly selects a fraction of the input structures.

        Args:
            structures: A list of ase.Atoms objects to sample from.

        Returns:
            A new list containing a random subset of the input structures.
        """
        if not structures:
            return []

        num_to_select = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_select)
