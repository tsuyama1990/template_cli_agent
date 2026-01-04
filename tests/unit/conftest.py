from typing import Any

import pytest

from mlip_autopipec.common.pydantic_models import (
    FullConfig,
)


@pytest.fixture
def valid_system_config_dict() -> dict[str, Any]:
    return {
        "elements": ["Cu", "Au"],
        "composition": {"Cu": 0.5, "Au": 0.5},
        "supercell_size": [3, 3, 3],
        "num_initial_structures": 10,
    }


@pytest.fixture
def valid_generation_config_dict() -> dict[str, Any]:
    return {
        "rattle_std_dev": 0.1,
        "volumetric_strain": 0.05,
        "min_atomic_distance": 2.0,
    }


@pytest.fixture
def valid_exploration_config_dict() -> dict[str, Any]:
    return {
        "md_calculator": "mace",
        "ensemble": "nvt",
        "temperature_k": 300.0,
        "time_step_fs": 1.0,
        "num_steps": 1000,
    }


@pytest.fixture
def valid_sampling_config_dict() -> dict[str, Any]:
    return {
        "method": "random",
        "num_samples": 100,
    }


@pytest.fixture
def valid_full_config_dict(
    valid_system_config_dict: dict[str, Any],
    valid_generation_config_dict: dict[str, Any],
    valid_exploration_config_dict: dict[str, Any],
    valid_sampling_config_dict: dict[str, Any],
) -> dict[str, Any]:
    return {
        "system": valid_system_config_dict,
        "generation": valid_generation_config_dict,
        "exploration": valid_exploration_config_dict,
        "sampling": valid_sampling_config_dict,
    }


@pytest.fixture
def mock_config(valid_full_config_dict: dict[str, Any]) -> FullConfig:
    """Provides a mock FullConfig object."""
    return FullConfig(**valid_full_config_dict)
