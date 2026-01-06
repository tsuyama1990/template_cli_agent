# -*- coding: utf-8 -*-
"""Unit tests for the Pydantic configuration models."""
from typing import Any

import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import (
    ExplorationConfig,
    FullConfig,
    SamplingConfig,
    SystemConfig,
)


def test_valid_system_config() -> None:
    """Test that a valid system configuration is parsed correctly."""
    config: dict[str, Any] = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": 10,
    }
    system_config = SystemConfig(**config)
    assert system_config.elements == ["Fe", "Pt"]
    assert system_config.composition == {"Fe": 0.5, "Pt": 0.5}
    assert system_config.lattice == "fcc"
    assert system_config.num_structures == 10


def test_invalid_composition_sum() -> None:
    """Test that a composition that does not sum to 1.0 raises an error."""
    config: dict[str, Any] = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.6, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": 10,
    }
    with pytest.raises(ValidationError, match="Composition must sum to 1.0"):
        SystemConfig(**config)


def test_mismatched_elements() -> None:
    """Test that mismatched elements and composition keys raise an error."""
    config: dict[str, Any] = {
        "elements": ["Fe", "Au"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": 10,
    }
    with pytest.raises(ValidationError, match="Elements and composition keys must match."):
        SystemConfig(**config)


def test_negative_num_structures() -> None:
    """Test that a negative number of structures raises an error."""
    config: dict[str, Any] = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": -1,
    }
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        SystemConfig(**config)


def test_valid_exploration_config() -> None:
    """Test that a valid exploration configuration is parsed correctly."""
    config = {"temperature": 300.0}
    exp_config = ExplorationConfig(**config)
    assert exp_config.temperature == 300.0


def test_negative_temperature() -> None:
    """Test that a negative temperature raises an error."""
    config = {"temperature": -100.0}
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        ExplorationConfig(**config)


def test_valid_sampling_config() -> None:
    """Test that a valid sampling configuration is parsed correctly."""
    config: dict[str, Any] = {"method": "random", "fraction": 0.8}
    sampling_config = SamplingConfig(**config)
    assert sampling_config.method == "random"
    assert sampling_config.fraction == 0.8


def test_invalid_sampling_fraction() -> None:
    """Test that a sampling fraction > 1.0 raises an error."""
    config: dict[str, Any] = {"method": "random", "fraction": 1.1}
    with pytest.raises(ValidationError, match="Input should be less than or equal to 1"):
        SamplingConfig(**config)


def test_invalid_sampling_method() -> None:
    """Test that an invalid sampling method raises an error."""
    config: dict[str, Any] = {"method": "invalid", "fraction": 0.8}
    with pytest.raises(ValidationError):
        SamplingConfig(**config)


def test_valid_full_config() -> None:
    """Test that a full, valid configuration is parsed correctly."""
    config: dict[str, Any] = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }
    full_config = FullConfig(**config)
    assert full_config.project_name == "test_project"
    assert full_config.system.lattice == "fcc"
    assert full_config.exploration.temperature == 300.0
    assert full_config.sampling.fraction == 0.8
