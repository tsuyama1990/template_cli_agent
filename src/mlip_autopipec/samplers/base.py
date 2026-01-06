from abc import ABC, abstractmethod

from ase import Atoms


class BaseSampler(ABC):
    """Abstract base class for all sampling algorithms."""

    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Selects a subset of structures from a trajectory.

        This method must be implemented by all concrete sampler classes.
        """
        raise NotImplementedError
