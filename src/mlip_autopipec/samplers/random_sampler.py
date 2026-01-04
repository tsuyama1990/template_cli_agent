import random

from ase import Atoms

from mlip_autopipec.common.pydantic_models import SamplingConfig

from .base import BaseSampler


class RandomSampler(BaseSampler):
    def __init__(self, config: SamplingConfig):
        self.config = config

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Randomly samples a specified number of frames from a trajectory.
        """
        if len(trajectory) < self.config.num_samples:
            raise ValueError(
                "Number of samples requested is greater than the trajectory length."
            )
        return random.sample(trajectory, self.config.num_samples)
