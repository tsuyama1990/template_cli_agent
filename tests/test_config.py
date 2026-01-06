import pytest
import yaml
from pydantic import ValidationError

from mlip_autopipec.config.models import FullConfig


def test_valid_config(tmp_path):
    """Test that a valid configuration file is loaded correctly."""
    config_content = """
project_name: test_project
system:
  elements: [Fe, Pt]
  composition: {Fe: 0.5, Pt: 0.5}
  lattice: fcc
  num_structures: 10
exploration:
  temperature: 300.0
sampling:
  method: random
  fraction: 0.5
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    with open(config_file, "r") as f:
        raw_config = yaml.safe_load(f)
    config = FullConfig.model_validate(raw_config)

    assert config.project_name == "test_project"
    assert config.system.elements == ["Fe", "Pt"]
    assert config.system.composition == {"Fe": 0.5, "Pt": 0.5}
    assert config.system.lattice == "fcc"
    assert config.system.num_structures == 10
    assert config.exploration.temperature == 300.0
    assert config.sampling.method == "random"
    assert config.sampling.fraction == 0.5


def test_invalid_composition_sum():
    """Test that a validation error is raised for invalid composition sum."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Pt"], "composition": {"Fe": 0.6, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_mismatched_elements():
    """Test that a validation error is raised for mismatched elements."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Au"], "composition": {"Fe": 0.5, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_negative_temperature():
    """Test that a validation error is raised for negative temperature."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Pt"], "composition": {"Fe": 0.5, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": -100},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)
