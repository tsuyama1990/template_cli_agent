import pytest
import os
import yaml
import numpy as np
import subprocess
from click.testing import CliRunner
from unittest.mock import MagicMock, patch

from mlip_autopipec.cli import main as cli_main
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.config import DFTResult
import logging

# A sample QE output file content
SAMPLE_QE_OUTPUT = """
     Final energy is     -450.0 eV
     Forces:
     atom    1   fx = 0.1   fy = 0.0   fz = 0.0
     atom    2   fx = -0.1  fy = 0.0   fz = 0.0
     atom    3   fx = 0.0   fy = 0.0   fz = 0.0

     Total stress:
     0.1 0.0 0.0
     0.0 0.1 0.0
     0.0 0.0 0.1
"""

@pytest.fixture
def mock_config_file(tmp_path):
    """Creates a temporary config file for testing."""
    config_data = {
        "database_path": str(tmp_path / "test.db"),
        "qe_command": "/path/to/mock/pw.x -in qe_calc.in > qe_calc.out",
        "dft_compute": {
            "pseudopotentials": {"H": "H.UPF", "O": "O.UPF"},
            "kpoints": [1, 1, 1],
            "ecutwfc": 50,
            "control": {"calculation": "scf"},
        },
        "mlip_training": {
            "model_type": "ACE",
            "r_cut": 4.5,
            "loss_weights": {"energy": 1.0, "forces": 1.0},
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)
    return config_path

@patch("mlip_autopipec.modules.labeling_engine.Espresso")
def test_end_to_end_workflow(mock_espresso_cls, mock_config_file, tmp_path):
    """
    Tests the full CLI workflow from db initialization to training.
    Mocks the external Quantum Espresso calculator.
    """
    runner = CliRunner()
    db_path = tmp_path / "test.db"

    # 1. Initialize the database
    result_init = runner.invoke(cli_main, ["init-db", "--db", str(db_path)])
    assert result_init.exit_code == 0
    assert "Added H2O molecule with ID: 1" in result_init.output

    # Configure the mock Espresso calculator
    mock_calc_instance = MagicMock()
    mock_espresso_cls.return_value = mock_calc_instance

    # Mock the results that would be read after the calculation
    mock_calc_instance.get_potential_energy.return_value = -450.0
    mock_calc_instance.get_forces.return_value = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0], [0.0, 0.0, 0.0]])
    mock_calc_instance.get_stress.return_value = np.eye(3) * 0.1

    # 2. Run the labeling command
    with patch('mlip_autopipec.modules.labeling_engine.logger') as mock_logger:
        result_label = runner.invoke(
            cli_main, ["label", "--id", "1", "--config", str(mock_config_file)],
            catch_exceptions=False # Show full traceback on error
        )
        assert result_label.exit_code == 0, result_label.output

        # Check that the logger was called with the success message
        mock_logger.info.assert_any_call("Successfully labeled structure ID: 1 with Energy: -450.0000 eV")

    # Verify that the database was updated correctly
    db_wrapper = AseDBWrapper(str(db_path))
    row = db_wrapper.get_row(1)
    assert row.key_value_pairs["state"] == "labeled"
    dft_result = DFTResult.model_validate_json(row.key_value_pairs["dft_result"])
    assert dft_result.energy == -450.0
    np.testing.assert_allclose(dft_result.forces[0, 0], 0.1)

    # 3. Run the training command
    result_train = runner.invoke(
        cli_main, ["train", "--config", str(mock_config_file)]
    )
    assert result_train.exit_code == 0, result_train.output
    assert "Training complete. Model saved to 'ace_model.json'" in result_train.output

    # Check that the dummy model file was created
    model_file_path = os.path.join(os.getcwd(), "ace_model.json")
    assert os.path.exists(model_file_path)
    os.remove(model_file_path) # Clean up
