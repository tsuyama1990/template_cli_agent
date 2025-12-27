"""Unit tests for the Labelling Engine."""
import pytest
from unittest.mock import patch, MagicMock
from ase.atoms import Atoms
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
import os

@pytest.fixture
def db_path(tmp_path):
    return os.path.join(tmp_path, "test.db")

@pytest.fixture
def mock_db(db_path):
    db = AseDB(db_path)
    db.write = MagicMock(return_value=1)
    return db

@pytest.fixture
def sample_atoms():
    return Atoms('H', positions=[[0, 0, 0]])

@pytest.fixture
def sample_config():
    return {
        "qe_command": "pw.x",
        "parameters": {
            "control": {"calculation": "'scf'"},
        },
        "pseudos": {"H": "H.upf"}
    }

@patch("subprocess.run")
def test_labelling_engine_success(mock_subprocess_run, mock_db, sample_atoms, sample_config):
    """Test the LabellingEngine for a successful QE run."""
    # Mock a successful QE run
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = "JOB DONE.! total energy = -1.0 Ry"
    mock_process.stderr = ""
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(
        qe_command=sample_config["qe_command"],
        parameters=sample_config["parameters"],
        pseudos=sample_config["pseudos"],
        db=mock_db
    )

    db_id = engine.execute(sample_atoms)

    # Verify that subprocess was called correctly
    mock_subprocess_run.assert_called_once()
    args, kwargs = mock_subprocess_run.call_args
    assert "pw.x -in" in args[0]

    # Verify that the db write method was called with a successful result
    mock_db.write.assert_called_once()
    called_atoms, called_result = mock_db.write.call_args[0]
    assert isinstance(called_atoms, Atoms)
    assert isinstance(called_result, DFTResult)
    assert called_result.was_successful
    assert db_id == 1

@patch("subprocess.run")
def test_labelling_engine_failure(mock_subprocess_run, mock_db, sample_atoms, sample_config):
    """Test the LabellingEngine for a failed QE run."""
    # Mock a failed QE run
    mock_process = MagicMock()
    mock_process.returncode = 1
    mock_process.stdout = "convergence has not been achieved"
    mock_process.stderr = "Error details"
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(
        qe_command=sample_config["qe_command"],
        parameters=sample_config["parameters"],
        pseudos=sample_config["pseudos"],
        db=mock_db
    )

    engine.execute(sample_atoms)

    # Verify that the db write method was called with a failed result
    mock_db.write.assert_called_once()
    _, called_result = mock_db.write.call_args[0]
    assert not called_result.was_successful
    assert "SCF did not converge" in called_result.error_message
