# tests/mlip_autopipec/integration/test_pipeline.py
"""Integration tests for the full CLI pipeline."""
import tempfile
import os
import yaml
import pytest
from typer.testing import CliRunner
from ase.db import connect
from mlip_autopipec.cli.main import app
import typer

runner = CliRunner()

def test_cli_pipeline_happy_path() -> None:
    """Test the full CLI pipeline with a valid configuration."""
    # 1. Create a temporary YAML config file
    config_data = {
        "system": {
            "elements": ["Cu"],
            "composition": {"Cu": 1.0},
            "supercell_size": [1, 1, 1],
        },
        "exploration": {"temperature_k": 300, "n_steps": 10},
        "sampling": {"method": "Random", "n_samples": 5},
    }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp_config:
        yaml.dump(config_data, tmp_config)
        config_path = tmp_config.name

    # Create a temporary path for the database
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp_db:
        db_path = tmp_db.name

    config_data["db_path"] = db_path
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)


    # 2. Run the CLI command
    result = runner.invoke(app, [config_path])

    # 3. Assertions
    assert result.exit_code == 0
    assert "Pipeline finished successfully!" in result.stdout

    # Verify the output database
    assert os.path.exists(db_path)
    with connect(db_path) as db:
        assert db.count() == 5  # Should match n_samples

    # Clean up the temporary files
    os.remove(config_path)
    os.remove(db_path)

def test_cli_handles_invalid_config() -> None:
    """Test that the CLI exits gracefully with an invalid config."""
    # 1. Create a config file with a missing required field ("elements")
    invalid_config_data = {
        "system": {
            "composition": {"Cu": 1.0},
            "supercell_size": [1, 1, 1],
        },
        "exploration": {"temperature_k": 300, "n_steps": 10},
        "sampling": {"n_samples": 5},
    }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp_config:
        yaml.dump(invalid_config_data, tmp_config)
        config_path = tmp_config.name

    # 2. Run the CLI command
    result = runner.invoke(app, [config_path])

    # 3. Assertions
    assert result.exit_code == 1
    assert "An error occurred" in result.stdout

    # Clean up
    os.remove(config_path)
