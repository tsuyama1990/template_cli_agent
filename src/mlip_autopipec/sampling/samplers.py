"""Implementation of sampling algorithms."""
import random
from typing import List

from ase.atoms import Atoms

from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """Selects a random subset of structures."""

    def __init__(self, config: SamplingConfig):
        """
        Initialize the sampler with sampling configuration.

        Args:
            config: The SamplingConfig Pydantic model.
        """
        self.config = config

    def sample(self, structures: List[Atoms]) -> List[Atoms]:
        """
        Randomly sample a fraction of the input structures.

        Args:
            structures: A list of ase.Atoms objects.

        Returns:
            A new list containing a random subset of the input structures.
        """
        num_to_select = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_select)
