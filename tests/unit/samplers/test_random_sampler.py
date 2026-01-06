"""Unit tests for the RandomSampler."""

import pytest
from mlip_autopipec.samplers.random import RandomSampler
from mlip_autopipec.common.pydantic_models import SamplingConfig
from ase import Atoms

@pytest.fixture
def sampling_config() -> SamplingConfig:
    """Provides a default SamplingConfig."""
    return SamplingConfig(method="random", n_samples=5)

def test_random_sampler_initialization(sampling_config: SamplingConfig) -> None:
    """Tests that the RandomSampler can be initialized."""
    sampler = RandomSampler(config=sampling_config)
    assert sampler is not None

def test_random_sampler_selects_correct_number(sampling_config: SamplingConfig) -> None:
    """
    Tests that the sampler selects the specified number of samples.
    """
    sampler = RandomSampler(config=sampling_config)
    # Create a dummy trajectory with more atoms than n_samples
    dummy_trajectory = [Atoms('H', positions=[(0, 0, 0)]) for _ in range(10)]

    selected = sampler.sample(dummy_trajectory)
    assert len(selected) == sampling_config.n_samples
