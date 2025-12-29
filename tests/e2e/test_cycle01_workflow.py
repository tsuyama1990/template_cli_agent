import os
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import run_cycle01_workflow


@pytest.fixture
def successful_qe_output():
    """Fixture to provide a sample successful QE output string."""
    with open("tests/unit/test_data/qe_outputs/h_atom_successful_run.out") as f:
        return f.read()

@patch('subprocess.run')
@patch('mlip_autopipec.modules.d_training_engine.TrainingEngine._save_model')
def test_e2e_cycle01_workflow_successful_run(
    mock_save_model, mock_subprocess_run, successful_qe_output
):
    """
    Tests the full end-to-end workflow for Cycle 01 using the CLI runner.
    Mocks the external QE call and the model saving to keep the test fast
    and dependency-free.
    """
    # Arrange
    runner = CliRunner()

    # Mock the QE execution
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = successful_qe_output
    mock_subprocess_run.return_value = mock_process

    # Mock the model saving to avoid actual training
    mock_save_model.return_value = "models/e2e_test_model.pt"

    config_path = "tests/e2e/test_config.yaml"
    structure_path = "tests/e2e/h_atom.xyz"
    db_path = "test_workflow.db"

    # Cleanup before test
    if os.path.exists(db_path):
        os.remove(db_path)

    # Act
    result = runner.invoke(
        run_cycle01_workflow,
        ['--config', config_path, '--structure', structure_path]
    )

    # Assert
    assert result.exit_code == 0
    assert "--- Starting MLIP-AutoPipe Cycle 01 Workflow ---" in result.output
    assert "Executing Labelling Engine..." in result.output
    assert "Labelling complete." in result.output
    assert "Executing Training Engine..." in result.output
    assert "Training complete. Model saved to: models/e2e_test_model.pt" in result.output
    assert "--- Workflow complete. ---" in result.output

    # Verify that the database was created
    assert os.path.exists(db_path)

    # Cleanup after test
    if os.path.exists(db_path):
        os.remove(db_path)
