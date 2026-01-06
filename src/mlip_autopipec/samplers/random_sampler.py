import random

from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.interfaces import BaseSampler


class RandomSampler(BaseSampler):
    """
    Selects a random subset of structures from a given trajectory.
    """

    def __init__(self, config: FullConfig):
        self.num_samples = config.sampling.num_samples

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Performs random sampling on the trajectory.

        Args:
            trajectory: A list of ASE Atoms objects from the simulation.

        Returns:
            A new list of ASE Atoms objects containing the sampled structures.
        """
        if not trajectory:
            return []

        if len(trajectory) <= self.num_samples:
            return trajectory  # Return all if trajectory is smaller than sample size

        return random.sample(trajectory, self.num_samples)
