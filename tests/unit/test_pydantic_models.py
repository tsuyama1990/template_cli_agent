from typing import Any

import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import (
    ExplorationConfig,
    FullConfig,
    GenerationConfig,
    SamplingConfig,
    SystemConfig,
)


def test_system_config_valid(valid_system_config_dict: dict[str, Any]) -> None:
    """Tests that a valid SystemConfig model can be created."""
    config = SystemConfig(**valid_system_config_dict)
    assert config.elements == ["Cu", "Au"]
    assert config.composition == {"Cu": 0.5, "Au": 0.5}


def test_system_config_composition_sum_invalid(
    valid_system_config_dict: dict[str, Any]
) -> None:
    """Tests that a validation error is raised if composition does not sum to 1.0."""
    valid_system_config_dict["composition"] = {"Cu": 0.5, "Au": 0.4}
    with pytest.raises(ValidationError, match=r"Composition values must sum to 1.0"):
        SystemConfig(**valid_system_config_dict)


def test_system_config_composition_keys_invalid(
    valid_system_config_dict: dict[str, Any]
) -> None:
    """Tests that a validation error is raised if composition keys are not in elements."""
    valid_system_config_dict["composition"] = {"Cu": 0.5, "Ag": 0.5}
    with pytest.raises(
        ValidationError, match=r"Composition keys must be a subset of elements"
    ):
        SystemConfig(**valid_system_config_dict)


def test_system_config_num_structures_invalid(
    valid_system_config_dict: dict[str, Any]
) -> None:
    """Tests that num_initial_structures must be positive."""
    valid_system_config_dict["num_initial_structures"] = 0
    with pytest.raises(ValidationError):
        SystemConfig(**valid_system_config_dict)


def test_generation_config_invalid(
    valid_generation_config_dict: dict[str, Any]
) -> None:
    """Tests that generation config raises errors for negative values."""
    valid_generation_config_dict["rattle_std_dev"] = -0.1
    with pytest.raises(ValidationError):
        GenerationConfig(**valid_generation_config_dict)

    valid_generation_config_dict["rattle_std_dev"] = 0.1
    valid_generation_config_dict["min_atomic_distance"] = 0
    with pytest.raises(ValidationError):
        GenerationConfig(**valid_generation_config_dict)


def test_exploration_config_invalid(
    valid_exploration_config_dict: dict[str, Any]
) -> None:
    """Tests that exploration config raises errors for non-positive values."""
    valid_exploration_config_dict["temperature_k"] = 0
    with pytest.raises(ValidationError):
        ExplorationConfig(**valid_exploration_config_dict)


def test_sampling_config_invalid(valid_sampling_config_dict: dict[str, Any]) -> None:
    """Tests that sampling config raises errors for non-positive values."""
    valid_sampling_config_dict["num_samples"] = 0
    with pytest.raises(ValidationError):
        SamplingConfig(**valid_sampling_config_dict)


def test_full_config_valid(valid_full_config_dict: dict[str, Any]) -> None:
    """Tests that a full, valid configuration can be loaded."""
    config = FullConfig(**valid_full_config_dict)
    assert config.system.elements == ["Cu", "Au"]
    assert config.sampling.num_samples == 100


def test_full_config_extra_fields_forbidden(
    valid_full_config_dict: dict[str, Any]
) -> None:
    """Tests that extra fields in the config raise a validation error."""
    valid_full_config_dict["system"]["unknown_field"] = "some_value"
    with pytest.raises(ValidationError, match=r"Extra inputs are not permitted"):
        FullConfig(**valid_full_config_dict)
