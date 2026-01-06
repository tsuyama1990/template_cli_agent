# src/mlip_autopipec/samplers/random_sampler.py
"""
Implements a simple random sampling strategy.
"""

import random

from ase import Atoms

from .base import BaseSampler


class RandomSampler(BaseSampler):
    """Selects a random subset of structures from a trajectory."""

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Performs random sampling on a trajectory.

        If the number of available structures is less than the number of
        requested samples, it returns all available structures.

        Args:
            trajectory: A list of ASE Atoms objects representing the trajectory.

        Returns:
            A new list of ASE Atoms objects containing the sampled structures.
        """
        n_samples = self.config.sampling.n_samples

        if len(trajectory) <= n_samples:
            return trajectory

        return random.sample(trajectory, n_samples)
