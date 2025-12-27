from pathlib import Path

import pytest

from mlip_autopipec.configs.config_loader import load_config


def test_load_config_success():
    """Tests that a valid configuration file is loaded correctly."""
    config_path = Path("tests/unit/test_config.yaml")
    config = load_config(config_path)
    assert config.training.model_path == "models/test_model.pt"
    assert config.training.epochs == 10
    assert config.training.delta_learn is True


def test_load_config_not_found():
    """Tests that a FileNotFoundError is raised for a missing config file."""
    with pytest.raises(FileNotFoundError):
        load_config(Path("non_existent_config.yaml"))


def test_load_config_invalid_schema():
    """Tests that a ValueError is raised for a config file with an invalid schema."""
    invalid_config_content = """
    training:
      model_path: "models/test_model.pt"
      # Missing 'epochs'
      delta_learn: true
    """
    invalid_config_path = Path("tests/unit/invalid_config.yaml")
    with open(invalid_config_path, "w") as f:
        f.write(invalid_config_content)

    with pytest.raises(ValueError):
        load_config(invalid_config_path)
    invalid_config_path.unlink()
