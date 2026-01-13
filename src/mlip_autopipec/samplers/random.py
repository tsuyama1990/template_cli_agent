import random

from ase import Atoms

from mlip_autopipec.samplers.base import BaseSampler


class RandomSampler(BaseSampler):
    """A sampler that randomly selects a subset of structures."""

    def __init__(self, num_samples: int = 3) -> None:
        self.num_samples = num_samples

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Randomly selects a specified number of structures from the trajectory.
        """
        if len(trajectory) <= self.num_samples:
            return trajectory
        return random.sample(trajectory, self.num_samples)
