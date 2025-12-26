import pytest
import yaml
from click.testing import CliRunner
from ase.io import write
from unittest.mock import patch

from mlip_autopipec.main_cycle01 import main

@pytest.fixture
def test_files(tmp_path):
    """Fixture to create necessary files for an E2E test."""
    # 1. Create a dummy structure file
    structure_path = tmp_path / "si.cif"
    from ase.build import bulk
    write(structure_path, bulk("Si", "diamond", a=5.43))

    # 2. Create a dummy config file
    config_data = {
        "qe_command": "pw.x",
        "db_path": str(tmp_path / "test.db"),
        "training": {
            "model_type": "mace",
            "learning_rate": 0.01,
            "num_epochs": 1,
            "r_cut": 4.0,
            "delta_learn": True,
            "baseline_potential": "lennard_jones",
        }
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    return config_path, structure_path

@patch("mlip_autopipec.modules.d_training_engine.run_train_cli")
@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
def test_e2e_cycle01_workflow_success(mock_subprocess_run, mock_train_cli, test_files):
    """
    Tests the full end-to-end workflow for a successful run.
    This satisfies UAT-C01-01.
    """
    config_path, structure_path = test_files

    # Mock the external process calls
    # 1. Mock Quantum Espresso to return a successful result
    mock_subprocess_run.return_value.returncode = 0
    mock_subprocess_run.return_value.stdout = """
    !    total energy              =     -150.0 Ry
    Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =    0.0 0.0 0.0
     atom    2 type  1   force =    0.0 0.0 0.0
    """
    # 2. The training CLI is already mocked

    runner = CliRunner()
    result = runner.invoke(main, ["--config", str(config_path), "--structure", str(structure_path)])

    # Verify the run was successful
    assert result.exit_code == 0
    assert "--- Starting Cycle 01 Workflow ---" in result.output
    assert "Executing Labelling Engine..." in result.output
    assert "Executing Training Engine..." in result.output
    assert "--- Workflow complete ---" in result.output
    assert "Trained model saved to:" in result.output

    # Verify that the mocked external calls were made
    mock_subprocess_run.assert_called_once()
    mock_train_cli.assert_called_once()
