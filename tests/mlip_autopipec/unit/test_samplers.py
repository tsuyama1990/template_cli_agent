# tests/mlip_autopipec/unit/test_samplers.py
"""Unit tests for sampling algorithms."""
from unittest.mock import Mock
from ase import Atoms
from mlip_autopipec.common.pydantic_models import SamplingConfig
from mlip_autopipec.samplers.random import RandomSampler


def test_random_sampler_returns_correct_number_of_samples() -> None:
    """Test that the RandomSampler returns the requested number of structures."""
    # 1. Configuration and mock trajectory
    config = SamplingConfig(method="Random", n_samples=5)
    trajectory = [Atoms("H") for _ in range(20)]  # A trajectory of 20 simple atoms

    # 2. Initialize sampler and run
    sampler = RandomSampler(config)
    sampled_structures = sampler.sample(trajectory)

    # 3. Assertions
    assert len(sampled_structures) == 5


def test_random_sampler_handles_small_trajectory() -> None:
    """Test the edge case where the trajectory is smaller than the number of samples."""
    # 1. Configuration and mock trajectory
    config = SamplingConfig(method="Random", n_samples=10)
    trajectory = [Atoms("He") for _ in range(5)]  # Trajectory is smaller than n_samples

    # 2. Initialize sampler and run
    sampler = RandomSampler(config)
    sampled_structures = sampler.sample(trajectory)

    # 3. Assertions
    # It should return the entire trajectory without errors
    assert len(sampled_structures) == 5
    assert sampled_structures == trajectory
