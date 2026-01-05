"""Implementation of sampling algorithms."""
import random

from ase import Atoms

from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """A simple random sampler."""

    def __init__(self, config: SamplingConfig):
        """Initialize the RandomSampler."""
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """Sample a random fraction of the structures."""
        num_to_sample = int(round(self.config.fraction * len(structures)))
        return random.sample(structures, num_to_sample)
