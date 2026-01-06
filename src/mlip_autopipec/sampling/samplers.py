from __future__ import annotations

import logging
import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import SamplingConfig

logger = logging.getLogger(__name__)


class RandomSampler:
    """A sampler that selects a random subset of structures."""

    def __init__(self, config: SamplingConfig):
        """
        Initialise the RandomSampler.

        Args:
            config: The sampling configuration.
        """
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Select a random subset of structures.

        Args:
            structures: The input structures.

        Returns:
            A randomly selected subset of the input structures.
        """
        num_to_select = int(len(structures) * self.config.fraction)
        sampled_structures = random.sample(structures, num_to_select)
        logger.info(f"Randomly sampled {len(sampled_structures)} structures.")
        return sampled_structures
