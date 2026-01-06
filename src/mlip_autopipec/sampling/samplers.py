from __future__ import annotations

import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """A simple sampler that randomly selects a fraction of structures."""

    def __init__(self, config: SamplingConfig):
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """Randomly selects a fraction of the input structures."""
        num_to_sample = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_sample)
