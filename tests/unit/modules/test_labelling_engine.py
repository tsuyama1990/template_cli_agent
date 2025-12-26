import pytest
from unittest.mock import MagicMock, patch
from ase.build import bulk
import subprocess
from pathlib import Path

from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.models import DFTResult

@pytest.fixture
def mock_db():
    """Fixture for a mocked AseDB instance."""
    db = MagicMock()
    db.write.return_value = 1  # Simulate returning a db_id
    return db

@pytest.fixture
def labelling_engine(mock_db):
    """Fixture for a LabellingEngine instance with a mocked database."""
    return LabellingEngine(
        qe_command="pw.x",
        db=mock_db,
        pseudos={"Si": "Si.upf"},
        kpts=(4, 4, 4),
        ecutwfc=80.0,
    )

@pytest.fixture
def sample_qe_output():
    """Loads the sample QE output file."""
    # We use the same sample output as the dft_utils test for consistency
    path = Path(__file__).parent.parent / "utils" / "qe_output.txt"
    with open(path, "r") as f:
        return f.read()


def test_labelling_engine_successful_run(labelling_engine, mock_db, sample_qe_output):
    """
    Tests the LabellingEngine's execute method for a successful QE run.
    """
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock the subprocess.run call to simulate a successful QE execution
    mock_process_result = MagicMock()
    mock_process_result.stdout = sample_qe_output

    with patch("subprocess.run", return_value=mock_process_result) as mock_run:
        db_id = labelling_engine.execute(atoms)

        # 1. Check that subprocess.run was called correctly
        mock_run.assert_called_once()
        assert "pw.x" in mock_run.call_args[0][0]

        # 2. Check that the database write method was called
        mock_db.write.assert_called_once()

        # 3. Check the content of what was written to the database
        written_result = mock_db.write.call_args[0][1]
        assert isinstance(written_result, DFTResult)
        assert written_result.was_successful is True
        assert written_result.total_energy_ev == pytest.approx(-11.24355702 * 13.605693122994)

        # 4. Check that the correct db_id was returned
        assert db_id == 1

def test_labelling_engine_subprocess_error(labelling_engine, mock_db):
    """
    Tests how the LabellingEngine handles a CalledProcessError from subprocess.
    """
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock subprocess.run to raise an error
    with patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "pw.x")) as mock_run:
        db_id = labelling_engine.execute(atoms)

        mock_run.assert_called_once()
        mock_db.write.assert_called_once()

        written_result = mock_db.write.call_args[0][1]
        assert written_result.was_successful is False
        assert "Subprocess execution failed" in written_result.error_message
        assert db_id == 1

def test_labelling_engine_qe_failure(labelling_engine, mock_db):
    """
    Tests how the LabellingEngine handles a failed (non-converged) QE run.
    """
    atoms = bulk("Si", "diamond", a=5.43)

    # Simulate a QE run that finished but did not converge
    mock_process_result = MagicMock()
    mock_process_result.stdout = " JOB DONE." # Missing convergence message

    with patch("subprocess.run", return_value=mock_process_result):
        labelling_engine.execute(atoms)

        written_result = mock_db.write.call_args[0][1]
        assert written_result.was_successful is False
        assert "SCF did not converge" in written_result.error_message
