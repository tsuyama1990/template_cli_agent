# tests/unit/samplers/test_random_sampler.py
import pytest
from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.samplers.random_sampler import RandomSampler


@pytest.fixture
def sampling_config() -> FullConfig:
    """Provides a valid config for random sampling."""
    config_dict = {
        "system": {
            "elements": ["H"],
            "composition": {"H": 1.0},
            "supercell_size": [1, 1, 1],
        },
        "exploration": {"temperature_k": 300, "pressure_gpa": 0},
        "sampling": {"method": "random", "n_samples": 5},
    }
    return FullConfig.model_validate(config_dict)


def test_random_sampler_selects_correct_number_of_samples(
    sampling_config: FullConfig,
) -> None:
    """Test that the sampler returns the requested number of structures."""
    # Create a dummy trajectory of 20 structures
    trajectory = [Atoms("H") for _ in range(20)]
    sampler = RandomSampler(config=sampling_config)

    sampled_structures = sampler.sample(trajectory)

    assert len(sampled_structures) == 5


def test_random_sampler_handles_fewer_structures_than_samples(
    sampling_config: FullConfig,
) -> None:
    """
    Test that the sampler returns all structures if the trajectory is smaller
    than the number of requested samples.
    """
    # Trajectory of 3 is less than the 5 requested samples
    trajectory = [Atoms("H") for _ in range(3)]
    sampler = RandomSampler(config=sampling_config)

    sampled_structures = sampler.sample(trajectory)

    # It should return all available structures without error
    assert len(sampled_structures) == 3
