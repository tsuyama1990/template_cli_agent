"""Unit tests for the configuration models."""
import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import FullConfig


def test_valid_config(tmp_path):
    """Test that a valid configuration is parsed correctly."""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        """
project_name: uat_test_fept
system:
  elements: ['Fe', 'Pt']
  composition: {'Fe': 0.5, 'Pt': 0.5}
  lattice: 'fcc'
  num_structures: 10
exploration:
  temperature: 300.0
sampling:
  method: 'random'
  fraction: 0.8
"""
    )
    import yaml

    with open(config_file) as f:
        config_dict = yaml.safe_load(f)

    config = FullConfig(**config_dict)
    assert config.project_name == "uat_test_fept"
    assert config.system.elements == ["Fe", "Pt"]
    assert config.system.composition == {"Fe": 0.5, "Pt": 0.5}
    assert config.system.lattice == "fcc"
    assert config.system.num_structures == 10
    assert config.exploration.temperature == 300.0
    assert config.sampling.method == "random"
    assert config.sampling.fraction == 0.8


def test_invalid_composition_sum():
    """Test that a validation error is raised for invalid composition sum."""
    with pytest.raises(ValidationError):
        FullConfig(
            project_name="test",
            system={
                "elements": ["Fe", "Pt"],
                "composition": {"Fe": 0.6, "Pt": 0.5},
                "lattice": "fcc",
                "num_structures": 10,
            },
            exploration={"temperature": 300.0},
            sampling={"method": "random", "fraction": 0.8},
        )


def test_mismatched_elements():
    """Test that a validation error is raised for mismatched elements."""
    with pytest.raises(ValidationError):
        FullConfig(
            project_name="test",
            system={
                "elements": ["Fe", "Pt"],
                "composition": {"Fe": 0.5, "Au": 0.5},
                "lattice": "fcc",
                "num_structures": 10,
            },
            exploration={"temperature": 300.0},
            sampling={"method": "random", "fraction": 0.8},
        )


def test_invalid_numeric_values():
    """Test that a validation error is raised for invalid numeric values."""
    with pytest.raises(ValidationError):
        FullConfig(
            project_name="test",
            system={
                "elements": ["Fe", "Pt"],
                "composition": {"Fe": 0.5, "Pt": 0.5},
                "lattice": "fcc",
                "num_structures": 10,
            },
            exploration={"temperature": -100.0},
            sampling={"method": "random", "fraction": 0.8},
        )

    with pytest.raises(ValidationError):
        FullConfig(
            project_name="test",
            system={
                "elements": ["Fe", "Pt"],
                "composition": {"Fe": 0.5, "Pt": 0.5},
                "lattice": "fcc",
                "num_structures": 10,
            },
            exploration={"temperature": 300.0},
            sampling={"method": "random", "fraction": 1.1},
        )


def test_incorrect_data_types():
    """Test that a validation error is raised for incorrect data types."""
    with pytest.raises(ValidationError):
        FullConfig(
            project_name="test",
            system={
                "elements": ["Fe", "Pt"],
                "composition": {"Fe": 0.5, "Pt": 0.5},
                "lattice": "fcc",
                "num_structures": "ten",
            },
            exploration={"temperature": 300.0},
            sampling={"method": "random", "fraction": 0.8},
        )
