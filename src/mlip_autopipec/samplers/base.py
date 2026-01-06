"""
Defines the base class for all sampling algorithms.
"""

from abc import ABC, abstractmethod

from ase import Atoms


class BaseSampler(ABC):
    """
    Abstract base class for all sampling algorithms.

    Each sampler must implement the `sample` method, which takes a trajectory
    of atomic structures and returns a selected subset.
    """

    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """
        Select a subset of structures from a trajectory.

        Args:
            trajectory: A list of ASE Atoms objects representing the trajectory.

        Returns:
            A list of ASE Atoms objects representing the selected samples.
        """
