"""Unit tests for the configuration models."""

import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import FullConfig


def test_valid_config() -> None:
    # Corresponds to SPEC.md, Section 3: Pydantic Schema Design
    # Verifies that a known-valid configuration object can be created.
    # Also part of UAT-C1-001: Successful End-to-End Run.
    """Test that a valid configuration is parsed correctly."""
    config_dict = {
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
    config = FullConfig(**config_dict)
    assert config.project_name == "test_project"
    assert config.system.lattice == "fcc"


def test_invalid_composition_sum() -> None:
    # Corresponds to SPEC.md, Section 3: Pydantic Schema Design (`@model_validator`)
    # Verifies that compositions not summing to 1.0 are rejected.
    # Also part of UAT-C1-003: Application Rejects a Configuration with Invalid Schema.
    """Test that a validation error is raised for invalid composition sum."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.6, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_mismatched_elements() -> None:
    # Corresponds to SPEC.md, Section 3: Pydantic Schema Design (`@model_validator`)
    # Verifies that a mismatch between `elements` and `composition` keys is rejected.
    """Test that a validation error is raised for mismatched elements."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Au"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_negative_temperature() -> None:
    # Corresponds to SPEC.md, Section 3: `ExplorationConfig`
    # Verifies the `gt=0` constraint on the temperature field.
    """Test that a validation error is raised for a negative temperature."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": -100.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_invalid_sampling_fraction() -> None:
    # Corresponds to SPEC.md, Section 3: `SamplingConfig`
    # Verifies the `gt=0` and `le=1` constraints on the fraction field.
    """Test that a validation error is raised for an invalid sampling fraction."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 1.1},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)
