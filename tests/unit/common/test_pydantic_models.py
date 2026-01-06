"""Unit tests for the Pydantic models in common.pydantic_models."""

import pytest
from pydantic import ValidationError
from typing import Dict, Any

from mlip_autopipec.common.pydantic_models import (
    SystemConfig,
    MDConfig,
    SamplingConfig,
    FullConfig,
)


def test_system_config_valid() -> None:
    """Tests that a valid SystemConfig model can be created."""
    config: Dict[str, Any] = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.25},
        "supercell_size": [3, 3, 3],
    }
    model = SystemConfig(**config)
    assert model.elements == ["Fe", "Pt"]
    assert model.supercell_size == [3, 3, 3]


def test_system_config_invalid_supercell() -> None:
    """Tests that SystemConfig raises an error for an invalid supercell size."""
    config: Dict[str, Any] = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.25},
        "supercell_size": [3, 3],  # Too short
    }
    with pytest.raises(ValidationError):
        SystemConfig(**config)


def test_md_config_valid() -> None:
    """Tests that a valid MDConfig model can be created."""
    config = {"temperature_k": 300.0, "pressure_gpa": 1.0}
    model = MDConfig(**config)
    assert model.temperature_k == 300.0


def test_md_config_invalid_temperature() -> None:
    """Tests that MDConfig raises an error for a non-positive temperature."""
    config = {"temperature_k": 0, "pressure_gpa": 1.0}
    with pytest.raises(ValidationError):
        MDConfig(**config)


def test_sampling_config_valid() -> None:
    """Tests that a valid SamplingConfig model can be created."""
    config: Dict[str, Any] = {"method": "random", "n_samples": 100}
    model = SamplingConfig(**config)
    assert model.method == "random"
    assert model.n_samples == 100


def test_full_config_valid() -> None:
    """Tests that a full, valid configuration can be parsed."""
    config: Dict[str, Any] = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 500.0, "pressure_gpa": 0.0},
        "sampling": {"method": "random", "n_samples": 50},
    }
    model = FullConfig(**config)
    assert model.system.elements == ["Si"]
    assert model.exploration.temperature_k == 500.0
    assert model.sampling.n_samples == 50


def test_full_config_extra_fields() -> None:
    """Tests that the FullConfig model forbids extra fields."""
    config: Dict[str, Any] = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 500.0, "pressure_gpa": 0.0},
        "sampling": {"method": "random", "n_samples": 50},
        "unexpected_field": "some_value",
    }
    with pytest.raises(ValidationError):
        FullConfig(**config)
