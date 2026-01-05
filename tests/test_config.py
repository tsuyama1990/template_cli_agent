"""Unit tests for the Pydantic configuration models."""

import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import (
    FullConfig,
)


def test_valid_config():
    """Test that a valid configuration is parsed correctly."""
    valid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    config = FullConfig(**valid_data)
    assert config.project_name == "test_project"
    assert config.system.lattice == "fcc"
    assert config.sampling.fraction == 0.8


def test_invalid_composition_sum():
    """Test that a validation error is raised for invalid composition sum."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.6, "Pt": 0.5},  # Sums to 1.1
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_mismatched_elements():
    """Test that a validation error is raised for mismatched elements."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Au": 0.5},  # Au not in elements
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_negative_temperature():
    """Test that a validation error is raised for negative temperature."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": -100.0},  # Invalid temperature
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_invalid_sampling_fraction():
    """Test that a validation error is raised for invalid sampling fraction."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 1.1},  # Invalid fraction
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_empty_elements_list():
    """Test that a validation error is raised for an empty elements list."""
    invalid_data = {
        "system": {
            "elements": [],
            "composition": {},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_zero_num_structures():
    """Test that a validation error is raised for zero num_structures."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 0,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)


def test_invalid_lattice_type():
    """Test that a validation error is raised for an invalid lattice type."""
    invalid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "invalid_lattice",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    with pytest.raises(ValidationError):
        FullConfig(**invalid_data)
