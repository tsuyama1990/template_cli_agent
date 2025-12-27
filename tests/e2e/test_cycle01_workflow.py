import pytest
import yaml
from click.testing import CliRunner
from unittest.mock import patch
import os

from mlip_autopipec.main_cycle01 import main

# A minimal, valid YAML config for the test
MOCK_CONFIG = {
    "db_path": "test_e2e.db",
    "qe_command": "mock_pw.x",
    "pseudos": {"Si": "Si.upf"},
    "dft_params": {
        "control": {"calculation": "scf"},
        "system": {"ecutwfc": 60.0},
        "electrons": {"mixing_beta": 0.7}
    },
    "training": {
        "model_type": "mace",
        "learning_rate": 0.001,
        "num_epochs": 1,
        "r_cut": 4.0,
        "delta_learn": True,
        "baseline_potential": "lj"
    }
}

# A minimal structure file (can be invalid CIF, just needs to exist)
MOCK_STRUCTURE_CONTENT = """
data_Si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
"""

@pytest.fixture
def test_environment(tmp_path):
    """Sets up a temporary directory with mock config and structure files."""
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(MOCK_CONFIG, f)

    structure_file = tmp_path / "structure.cif"
    with open(structure_file, "w") as f:
        f.write(MOCK_STRUCTURE_CONTENT)

    # Also create a mock pseudo file
    pseudo_file = tmp_path / "Si.upf"
    with open(pseudo_file, "w") as f:
        f.write("mock pseudo")

    # Change to the temp directory so that QE finds the pseudo
    os.chdir(tmp_path)

    yield config_file, structure_file

    # Cleanup
    os.remove("test_e2e.db") if os.path.exists("test_e2e.db") else None
    model_file = "models/trained_model.pt"
    if os.path.exists(model_file):
         os.remove(model_file)
         os.rmdir("models")


@patch('mlip_autopipec.modules.d_training_engine.TrainingEngine.execute')
@patch('mlip_autopipec.modules.c_labelling_engine.LabellingEngine.execute')
def test_e2e_workflow_successful(mock_labeller_execute, mock_trainer_execute, test_environment):
    """
    Tests the full end-to-end workflow from the CLI, mocking the engine executions.
    """
    config_file, structure_file = test_environment

    # Configure mocks to return successful results
    mock_labeller_execute.return_value = 1 # Mock database ID
    mock_trainer_execute.return_value = "models/trained_model.pt"

    # Mock the AseDB.get method to simulate a successful labelling
    with patch('mlip_autopipec.data.database.AseDB.get') as mock_db_get:
        from ase import Atoms
        mock_db_get.return_value = (Atoms('H'), {'was_successful': True})

        runner = CliRunner()
        result = runner.invoke(main, ['--config', str(config_file), '--structure', str(structure_file)])

    assert result.exit_code == 0
    assert "--- Starting Cycle 01 Workflow ---" in result.output
    assert "Labelling complete" in result.output
    assert "Workflow complete" in result.output

    mock_labeller_execute.assert_called_once()
    mock_trainer_execute.assert_called_once()


@patch('mlip_autopipec.modules.c_labelling_engine.LabellingEngine.execute')
def test_e2e_workflow_labelling_fails(mock_labeller_execute, test_environment):
    """
    Tests that the workflow aborts correctly if the labelling engine fails.
    """
    config_file, structure_file = test_environment
    mock_labeller_execute.return_value = 1

    with patch('mlip_autopipec.data.database.AseDB.get') as mock_db_get:
        from ase import Atoms
        mock_db_get.return_value = (Atoms('H'), {'was_successful': False, 'error_message': 'Test Failure'})

        runner = CliRunner()
        result = runner.invoke(main, ['--config', str(config_file), '--structure', str(structure_file)])

    assert result.exit_code == 0
    assert "ERROR: Labelling failed" in result.output
    assert "Aborting workflow" in result.output
    assert "Running Training Engine" not in result.output
