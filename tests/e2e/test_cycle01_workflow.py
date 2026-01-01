from unittest.mock import patch, MagicMock
import os
import yaml

import pytest
from ase import Atoms
from click.testing import CliRunner

# Under test
from mlip_autopipec.cli import mlip_pipe
from mlip_autopipec.workflow import WorkflowOrchestrator

# Sample data
SAMPLE_CONFIG = {
    "system": {"elements": ["Si"]},
    "dft_compute": {
        "code": "quantum_espresso",
        "command": "mpirun -np 1 pw.x",
        "pseudopotentials": "SSSP",
        "ecutwfc": 50.0,
        "ecutrho": 200.0,
        "kpoints_density": 3.0,
        "smearing": "mv",
        "degauss": 0.01,
    },
    "mlip_training": {
        "model_type": "ace",
        "r_cut": 5.5,
        "delta_learning": False,
        "base_potential": "lj_auto",
        "loss_weights": {"energy": 1.0, "force": 10.0},
    },
}

PARSED_QE_DATA = {
    "energy": -100.0,
    "forces": [[0.0] * 3],
    "stress": [[0.0] * 3] * 3,
}
SAMPLE_QE_OUTPUT_STR = "Fake QE output"

@pytest.fixture
def cli_runner(tmp_path):
    """Fixture to create a temporary directory and run the CLI inside it."""
    # Create temp files needed for the run
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(SAMPLE_CONFIG, f)

    db_path = tmp_path / "test.db"
    input_path = tmp_path / "inputs.xyz"

    # Create a dummy input file with one atom
    atoms = Atoms("Si", positions=[[0, 0, 0]])
    atoms.write(input_path)

    runner = CliRunner()

    # Yield the runner and paths to the test
    yield runner, config_path, db_path, input_path

    # Teardown (handled by tmp_path)


# We patch the dependencies of the entire workflow
@patch("mlip_autopipec.workflow.AseDBWrapper")
@patch("mlip_autopipec.workflow.LabelingEngine")
@patch("mlip_autopipec.workflow.TrainingEngine")
def test_workflow_orchestrator_integration(
    mock_training_engine, mock_labeling_engine, mock_db_wrapper, cli_runner
):
    """
    Test the WorkflowOrchestrator's ability to instantiate and run the engines.
    This is an integration test for the orchestrator, not a full E2E test.
    """
    runner, config_path, db_path, input_path = cli_runner

    # Instantiate mocks for the engines and db
    mock_db_instance = MagicMock()
    mock_labeling_instance = MagicMock()
    mock_training_instance = MagicMock()

    # The constructor of the orchestrator will receive these mocked classes
    # and should instantiate them
    mock_db_wrapper.return_value = mock_db_instance
    mock_labeling_engine.return_value = mock_labeling_instance
    mock_training_engine.return_value = mock_training_instance

    # Create the orchestrator with the path to the config file
    orchestrator = WorkflowOrchestrator(
        config_path=str(config_path),
        db_path=str(db_path),
        input_file_path=str(input_path)
    )

    # Run the main workflow
    orchestrator.run_workflow()

    # --- Assertions ---
    # 1. DB was connected to and atoms were added
    mock_db_instance.connect.assert_called_once_with(str(db_path))
    mock_db_instance.add_atoms.assert_called_once()

    # 2. LabelingEngine was instantiated with the correct config and db
    mock_labeling_engine.assert_called_once()
    call_args = mock_labeling_engine.call_args[1]
    assert call_args['config'].code == "quantum_espresso"
    assert call_args['db_wrapper'] is mock_db_instance

    # 3. LabelingEngine's run method was called
    mock_labeling_instance.run.assert_called_once()

    # 4. TrainingEngine was instantiated correctly
    mock_training_engine.assert_called_once()
    call_args = mock_training_engine.call_args[1]
    assert call_args['config'].model_type == "ace"
    assert call_args['db_wrapper'] is mock_db_instance

    # 5. TrainingEngine's run method was called
    mock_training_instance.run.assert_called_once()


@patch("mlip_autopipec.cli.WorkflowOrchestrator")
def test_cli_e2e(mock_orchestrator, cli_runner):
    """
    Test the CLI entrypoint. This is a higher-level E2E test.
    We mock the entire orchestrator to isolate the test to the CLI logic.
    """
    runner, config_path, db_path, input_path = cli_runner
    mock_orchestrator_instance = MagicMock()
    mock_orchestrator.return_value = mock_orchestrator_instance

    result = runner.invoke(mlip_pipe, [
        "--config-file", str(config_path),
        "--database-file", str(db_path),
        "--input-file", str(input_path),
    ])

    # Assert the CLI exited successfully
    assert result.exit_code == 0
    assert "Workflow initialized" in result.output

    # Assert the orchestrator was created and run
    mock_orchestrator.assert_called_once_with(
        config_path=str(config_path),
        db_path=str(db_path),
        input_file_path=str(input_path),
    )
    mock_orchestrator_instance.run_workflow.assert_called_once()
