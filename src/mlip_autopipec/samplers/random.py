# src/mlip_autopipec/samplers/random.py
"""Concrete implementation of a random sampling algorithm."""

import random

from ase import Atoms

from mlip_autopipec.common.pydantic_models import SamplingConfig
from mlip_autopipec.samplers.base import BaseSampler


class RandomSampler(BaseSampler):
    """Selects a random subset of structures from a trajectory."""

    def __init__(self, config: SamplingConfig) -> None:
        self.config = config

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Performs random sampling on the provided trajectory.

        Args:
            trajectory: A list of ASE Atoms objects.

        Returns:
            A randomly selected subset of the trajectory.
        """
        if len(trajectory) <= self.config.n_samples:
            return trajectory

        return random.sample(trajectory, self.config.n_samples)
