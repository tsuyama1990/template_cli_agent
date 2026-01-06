import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import (
    FullConfig,
)


def test_valid_config():
    """Tests that a valid configuration is parsed correctly."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    config = FullConfig(**config_dict)
    assert config.project_name == "test_project"
    assert config.system.lattice == "fcc"
    assert config.exploration.temperature == 300.0


def test_invalid_composition_sum():
    """Tests that a validation error is raised for invalid composition sum."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.6, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_mismatched_elements():
    """Tests that a validation error is raised for mismatched elements."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Au"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_invalid_temperature():
    """Tests that a validation error is raised for a negative temperature."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": -100.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)


def test_invalid_sampling_fraction():
    """Tests that a validation error is raised for a sampling fraction > 1.0."""
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


def test_incorrect_data_type():
    """Tests that a validation error is raised for an incorrect data type."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": "ten",
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config_dict)
