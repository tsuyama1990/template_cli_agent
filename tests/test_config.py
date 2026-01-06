from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from mlip_autopipec.config.models import FullConfig


def test_valid_config(tmp_path: Path) -> None:
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

    with config_file.open("r") as f:
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


def test_invalid_composition_sum() -> None:
    """Test that a validation error is raised for invalid composition sum."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Pt"], "composition": {"Fe": 0.6, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_mismatched_elements() -> None:
    """Test that a validation error is raised for mismatched elements."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Au"], "composition": {"Fe": 0.5, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_negative_temperature() -> None:
    """Test that a validation error is raised for negative temperature."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe", "Pt"], "composition": {"Fe": 0.5, "Pt": 0.5}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": -100},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_empty_elements_list() -> None:
    """Test that a validation error is raised for an empty elements list."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": [], "composition": {}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_zero_num_structures() -> None:
    """Test that a validation error is raised for zero num_structures."""
    raw_config = {
        "project_name": "test",
        "system": {"elements": ["Fe"], "composition": {"Fe": 1.0}, "lattice": "fcc", "num_structures": 0},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config)


def test_invalid_sampling_fraction() -> None:
    """Test validation error for sampling fraction outside the (0, 1] range."""
    # Case 1: fraction > 1
    raw_config_over = {
        "project_name": "test",
        "system": {"elements": ["Fe"], "composition": {"Fe": 1.0}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 1.1},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config_over)

    # Case 2: fraction == 0
    raw_config_zero = {
        "project_name": "test",
        "system": {"elements": ["Fe"], "composition": {"Fe": 1.0}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300},
        "sampling": {"method": "random", "fraction": 0.0},
    }
    with pytest.raises(ValidationError):
        FullConfig.model_validate(raw_config_zero)
