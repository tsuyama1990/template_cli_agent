import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import yaml
import subprocess
import os

from mlip_autopipec.main_cycle01 import run_workflow
from mlip_autopipec.data.database import AseDB

@pytest.fixture
def test_data_path():
    """Returns the path to the E2E test data directory."""
    return Path(__file__).parent

@pytest.fixture
def config(test_data_path, tmp_path):
    """Loads the test config and updates the database path to a temporary directory."""
    config_path = test_data_path / "test_config.yaml"
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)

    # Ensure the database is created in a temporary directory for isolation
    cfg['database']['path'] = str(tmp_path / "test_e2e.db")

    return cfg

@pytest.fixture
def sample_qe_output(test_data_path):
    """Loads the sample QE output file from the utils test data."""
    path = test_data_path.parent / "unit/utils/qe_output.txt"
    with open(path, "r") as f:
        return f.read()

def test_e2e_successful_workflow(config, test_data_path, sample_qe_output, tmp_path, monkeypatch):
    """
    Tests the full end-to-end workflow for a successful run (UAT-C01-01 & UAT-C01-02).
    """
    structure_path = str(test_data_path / "Si_bulk.cif")

    mock_process_result = MagicMock()
    mock_process_result.stdout = sample_qe_output

    # Use monkeypatch to safely change the working directory for the test
    monkeypatch.chdir(tmp_path)

    with patch("subprocess.run", return_value=mock_process_result) as mock_run:
        run_workflow(config, structure_path)

        # 1. Verify that QE was called
        mock_run.assert_called_once()
        assert "pw.x" in mock_run.call_args.args[0]

        # 2. Verify the database contents (UAT-C01-02)
        db = AseDB(config['database']['path'])
        row = db.get(1)
        assert row.key_value_pairs["was_successful"] is True
        assert row.key_value_pairs["total_energy_ev"] < 0

        dft_result = row.data
        forces = dft_result['forces']
        assert max(abs(f) for force_vec in forces for f in force_vec) < 1.0

        # 3. Verify that the model file was created (UAT-C01-01)
        model_path = Path(f"models/{config['training']['model_type']}_model.pt")
        assert model_path.is_file()
        assert model_path.stat().st_size > 0


def test_e2e_failed_dft_workflow(config, test_data_path, tmp_path, monkeypatch):
    """
    Tests the workflow's handling of a failed DFT calculation (UAT-C01-04).
    """
    structure_path = str(test_data_path / "Si_bulk.cif")

    mock_process_result = MagicMock()
    mock_process_result.stdout = "JOB DONE." # No convergence message

    monkeypatch.chdir(tmp_path)

    with patch("subprocess.run", return_value=mock_process_result):
        run_workflow(config, structure_path)

        # 1. Verify database shows failure
        db = AseDB(config['database']['path'])
        row = db.get(1)
        assert row.key_value_pairs["was_successful"] is False

        # 2. Verify that NO model file was created
        model_path = Path(f"models/{config['training']['model_type']}_model.pt")
        assert not model_path.is_file()
