"""A concrete implementation for random sampling."""

import random
from typing import List
from ase import Atoms
import logging

from mlip_autopipec.samplers.base import BaseSampler
from mlip_autopipec.common.pydantic_models import SamplingConfig

# Configure logger for this module
logger = logging.getLogger(__name__)

class RandomSampler(BaseSampler):
    """
    Selects a random subset of structures from a trajectory.
    """

    def __init__(self, config: SamplingConfig) -> None:
        """
        Initializes the RandomSampler.

        Args:
            config: The sampling configuration Pydantic model.
        """
        self.config = config
        if self.config.method != "random":
            msg = "RandomSampler only supports the 'random' method."
            raise ValueError(msg)
        logger.info(f"Initialized RandomSampler with config: {config}")

    def sample(self, trajectory: List[Atoms]) -> List[Atoms]:
        """
        Performs random sampling on the given trajectory.

        If the number of requested samples (`n_samples`) is greater than or
        equal to the number of structures in the trajectory, the entire
        trajectory is returned.

        Args:
            trajectory: A list of ASE Atoms objects representing the trajectory.

        Returns:
            A list of ASE Atoms objects representing the randomly selected samples.
        """
        n_structures = len(trajectory)
        n_samples = self.config.n_samples
        logger.info(
            f"Starting random sampling on a trajectory of {n_structures} structures to select {n_samples} samples."
        )

        if n_samples >= n_structures:
            logger.warning(
                f"Requested {n_samples} samples, but trajectory only has "
                f"{n_structures} structures. Returning the entire trajectory."
            )
            return trajectory

        try:
            sampled_indices = random.sample(range(n_structures), n_samples) # nosec B311
        except ValueError as e:
            logger.exception("Error during random sampling:")
            # This can happen if the population is smaller than the sample size,
            # though the check above should prevent it.
            msg = "Error during random sampling"
            raise ValueError(msg) from e
        else:
            selected_structures = [trajectory[i] for i in sorted(sampled_indices)]
            logger.info(f"Successfully sampled {len(selected_structures)} structures.")
            return selected_structures
