from __future__ import annotations

import random
from typing import TYPE_CHECKING, Protocol

import numpy as np
from dscribe.descriptors import SOAP
from scipy.spatial.distance import pdist, squareform

if TYPE_CHECKING:
    from ase import Atoms


class Sampler(Protocol):
    """Protocol for a sampler class."""

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """Select a subset of structures from a list."""
        ...


class RandomSampler:
    """Selects a random subset of structures."""

    def __init__(self, fraction: float = 0.5):
        self.fraction = fraction

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Randomly select a fraction of the provided structures.

        Args:
            structures: A list of `ase.Atoms` objects.

        Returns:
            A randomly selected sub-list of `ase.Atoms` objects.
        """
        if self.fraction >= 1.0:
            return structures
        num_to_select = int(len(structures) * self.fraction)
        return random.sample(structures, num_to_select)


class FPSSampler:
    """
    Selects a diverse subset of structures using Farthest Point Sampling (FPS).

    FPS is an iterative algorithm that selects points that are maximally distant
    from the set of already selected points, ensuring a diverse and representative
    subset.
    """

    def __init__(self, num_samples: int):
        self.num_samples = num_samples
        # Setup SOAP descriptor generator. These parameters are a reasonable
        # starting point for many atomic systems.
        self.soap = SOAP(
            species=["H", "C", "N", "O", "Cu", "Au", "Si", "Pt", "Fe"],
            periodic=True,
            r_cut=4.0,
            n_max=6,
            l_max=5,
            average="inner",
            sparse=False,
        )

    def _get_feature_vectors(self, structures: list[Atoms]) -> np.ndarray:
        """Calculate SOAP feature vectors for all structures."""
        # Note: This can be parallelized for better performance on large datasets.
        return self.soap.create(structures, n_jobs=1)

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Select a subset of structures using Farthest Point Sampling.

        Args:
            structures: A list of `ase.Atoms` objects to sample from.

        Returns:
            A diverse sub-list of `ase.Atoms` objects.
        """
        if len(structures) <= self.num_samples:
            return structures

        features = self._get_feature_vectors(structures)
        n_structures = len(structures)

        # Calculate the pairwise distance matrix
        distance_matrix = squareform(pdist(features, "euclidean"))

        # Start with a random structure
        first_index = random.randint(0, n_structures - 1)
        selected_indices = {first_index}

        # Initialize distances from each point to the selected set
        min_distances = distance_matrix[:, first_index]

        while len(selected_indices) < self.num_samples:
            # Find the point that is farthest from the current selection
            next_index = np.argmax(min_distances)
            selected_indices.add(int(next_index))

            # Update the minimum distances
            min_distances = np.minimum(min_distances, distance_matrix[:, next_index])

        return [structures[i] for i in sorted(list(selected_indices))]
