from unittest.mock import MagicMock, patch

import pytest
import yaml
from ase.build import bulk
from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import main


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture
def setup_test_environment(tmp_path):
    """Creates a temporary test environment with config and structure files."""
    # Create a dummy structure file
    atoms = bulk("Si", "diamond", a=5.43)
    structure_file = tmp_path / "Si_bulk.xyz"
    atoms.write(structure_file)

    # Create a dummy config file
    config_data = {
        "quantum_espresso": {"command": "pw.x -in QE.in > QE.out"},
        "training": {
            "model_type": "mace",
            "learning_rate": 0.01,
            "num_epochs": 5,
            "r_cut": 4.0,
            "delta_learn": False,
            "baseline_potential": "zbl",
        },
    }
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    return config_file, structure_file


# Mock all the heavy external libraries and processes
@patch("mlip_autopipec.orchestrator_cycle01.LabellingEngine.execute")
@patch("mlip_autopipec.orchestrator_cycle01.TrainingEngine.execute")
@patch("mlip_autopipec.orchestrator_cycle01.AseDB")
def test_e2e_workflow_success(
    mock_db_class, mock_trainer_execute, mock_labeller_execute, runner, setup_test_environment
):
    """
    Tests the full end-to-end workflow from the CLI on a successful run.
    """
    config_file, structure_file = setup_test_environment

    # Mock behavior
    mock_db_instance = MagicMock()
    mock_db_instance.get.return_value = (MagicMock(), {"was_successful": True})
    mock_db_class.return_value = mock_db_instance

    mock_labeller_execute.return_value = 1
    mock_trainer_execute.return_value = "models/test_model.pt"

    result = runner.invoke(main, ["--config", str(config_file), "--structure", str(structure_file)])

    assert result.exit_code == 0
    assert "--- Workflow Finished Successfully ---" in result.output
    assert "Model saved to: models/test_model.pt" in result.output

    mock_labeller_execute.assert_called_once()
    mock_trainer_execute.assert_called_once_with(ids=[1])


@patch("mlip_autopipec.orchestrator_cycle01.LabellingEngine.execute")
@patch("mlip_autopipec.orchestrator_cycle01.TrainingEngine.execute")
@patch("mlip_autopipec.orchestrator_cycle01.AseDB")
def test_e2e_workflow_labelling_fails(
    mock_db_class, mock_trainer_execute, mock_labeller_execute, runner, setup_test_environment
):
    """
    Tests that the workflow correctly halts if the labelling step fails.
    """
    config_file, structure_file = setup_test_environment

    # Mock behavior for failed labelling
    mock_db_instance = MagicMock()
    mock_db_instance.get.return_value = (
        MagicMock(),
        {"was_successful": False, "error_message": "SCF failed"},
    )
    mock_db_class.return_value = mock_db_instance

    mock_labeller_execute.return_value = 1

    result = runner.invoke(main, ["--config", str(config_file), "--structure", str(structure_file)])

    assert result.exit_code == 0
    assert "--- Workflow Finished (with errors) ---" in result.output
    assert "Skipping training because labelling failed" in result.output

    mock_labeller_execute.assert_called_once()
    # The trainer should NOT have been called
    mock_trainer_execute.assert_not_called()
