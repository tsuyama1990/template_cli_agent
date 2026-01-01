import json
import unittest.mock
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from mlip_autopipec.data.models import FullConfig
from mlip_autopipec.services.config_expander import ConfigExpander


@pytest.fixture
def mock_sssp_defaults() -> dict:
    """Provides a mock dictionary of SSSP default values."""
    return {
        "Fe": {"filename": "Fe.upf", "ecutwfc": 80.0, "ecutrho": 640.0},
        "Pt": {"filename": "Pt.upf", "ecutwfc": 90.0, "ecutrho": 720.0},
    }


def test_expand_config_success(tmp_path: Path, mock_sssp_defaults: dict):
    """
    Tests that a valid minimal config is expanded correctly into a FullConfig.
    """
    # 1. Arrange: Create minimal input and mock defaults file
    input_yaml_path = tmp_path / "input.yaml"
    input_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": "FePt",
        }
    }
    with open(input_yaml_path, "w") as f:
        yaml.dump(input_data, f)

    defaults_path = tmp_path / "sssp_defaults.json"
    with open(defaults_path, "w") as f:
        json.dump(mock_sssp_defaults, f)

    # 2. Act: Run the expander
    expander = ConfigExpander(defaults_path=defaults_path)
    full_config = expander.expand_config(input_path=input_yaml_path)

    # 3. Assert: Check that the heuristics were applied correctly
    assert isinstance(full_config, FullConfig)

    # Assert system inference
    assert full_config.system.structure_type == "alloy"
    assert full_config.system.composition == {"Fe": 1, "Pt": 1}

    # Assert DFT parameter selection (should pick the max of Fe and Pt)
    assert full_config.dft_compute.ecutwfc == 90.0
    assert full_config.dft_compute.ecutrho == 720.0
    assert full_config.dft_compute.pseudopotentials == {
        "Fe": "Fe.upf",
        "Pt": "Pt.upf",
    }

    # Assert default values were populated
    assert full_config.generation.generation_strategy == "sqs"
    assert full_config.mlip_training.model_type == "mace"
    assert full_config.mlip_training.delta_learning is False


def test_expand_config_invalid_input(tmp_path: Path):
    """
    Tests that the ConfigExpander raises a ValidationError for invalid input.
    """
    # 1. Arrange: Create an invalid minimal input (e.g., missing 'elements')
    input_yaml_path = tmp_path / "input.yaml"
    input_data = {
        "system": {
            "composition": "FePt",
            # 'elements' field is missing
        }
    }
    with open(input_yaml_path, "w") as f:
        yaml.dump(input_data, f)

    defaults_path = tmp_path / "sssp_defaults.json"
    with open(defaults_path, "w") as f:
        json.dump({}, f)

    # 2. Act & Assert: Expect a Pydantic validation error
    expander = ConfigExpander(defaults_path=defaults_path)
    with pytest.raises(ValidationError):
        expander.expand_config(input_path=input_yaml_path)
