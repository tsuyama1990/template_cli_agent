"""End-to-end test for the Cycle 01 workflow."""
import pytest
import os
import yaml
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
from mlip_autopipec.main import main

@pytest.fixture
def test_files(tmp_path):
    """Fixture to create necessary files for an e2e test."""
    config = {
        "qe_command": "pw.x",
        "database_path": os.path.join(tmp_path, "test.db"),
        "pseudopotentials": {"H": "H.UPF"},
        "dft_parameters": {
            "control": {"calculation": "'scf'"},
        },
        "training": {
            "model_type": "dummy",
            "learning_rate": 0.01,
            "num_epochs": 1,
            "r_cut": 4.0,
            "delta_learn": False,
            "baseline_potential": "none"
        }
    }
    config_path = os.path.join(tmp_path, "config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    structure_content = "2\nH2\nH 0 0 0\nH 0 0 0.74"
    structure_path = os.path.join(tmp_path, "h2.xyz")
    with open(structure_path, "w") as f:
        f.write(structure_content)

    return config_path, structure_path

@patch("subprocess.run")
def test_e2e_workflow_dummy_trainer(mock_subprocess_run, test_files):
    """Test the full e2e workflow with the dummy trainer."""
    config_path, structure_path = test_files

    # Mock a successful QE run
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = "JOB DONE."
    mock_subprocess_run.return_value = mock_process

    runner = CliRunner()
    result = runner.invoke(main, ['--config', config_path, '--structure', structure_path])

    assert result.exit_code == 0
    assert "Workflow Finished Successfully" in result.output

    # Check that the dummy model file was created
    model_path = os.path.join(os.getcwd(), "models", "dummy_model.pt")
    assert os.path.exists(model_path)

    # Clean up
    os.remove(model_path)
    os.rmdir("models")
