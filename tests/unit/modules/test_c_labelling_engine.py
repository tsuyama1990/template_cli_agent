"""Unit tests for the LabellingEngine module."""
import pytest
import subprocess
from unittest.mock import MagicMock, patch
from pathlib import Path
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine

# Re-use the sample QE outputs from the dft_utils test
SAMPLE_DATA_DIR = Path(__file__).parent.parent / "utils/sample_qe_outputs"

@pytest.fixture
def mock_db(tmp_path: Path) -> AseDB:
    """Provides a mock AseDB instance with a temporary db file."""
    db_path = tmp_path / "test.db"
    return AseDB(db_path)

@pytest.fixture
def labelling_engine(mock_db: AseDB) -> LabellingEngine:
    """Provides a configured LabellingEngine instance."""
    return LabellingEngine(
        db=mock_db,
        qe_command="pw.x",
        pseudo_dir="/path/to/pseudos",
        pseudos={"Si": "Si.upf"},
        kpts=(4, 4, 4),
        ecutwfc=60.0,
    )

@patch("subprocess.run")
def test_labelling_engine_execute_success(mock_subprocess_run, labelling_engine: LabellingEngine, mock_db: AseDB):
    """
    Tests a successful run of the LabellingEngine's execute method.
    """
    # Arrange
    atoms = bulk("Si", "diamond", a=5.43)
    success_output = (SAMPLE_DATA_DIR / "success.out").read_text()

    # Configure the mock to return a successful process
    mock_subprocess_run.return_value = subprocess.CompletedProcess(
        args=["pw.x", "-in", "qe.in"],
        returncode=0,
        stdout=success_output,
        stderr="",
    )

    # Act
    db_id = labelling_engine.execute(atoms)

    # Assert
    mock_subprocess_run.assert_called_once() # Check that QE was called

    # Verify the data was written to the DB correctly
    retrieved_atoms, result = mock_db.get(db_id)
    assert result.was_successful is True
    assert result.total_energy_ev is not None
    assert len(retrieved_atoms) == 2

@patch("subprocess.run")
def test_labelling_engine_execute_qe_failure(mock_subprocess_run, labelling_engine: LabellingEngine, mock_db: AseDB):
    """
    Tests the LabellingEngine's handling of a failed QE run (non-zero exit code).
    """
    # Arrange
    atoms = bulk("Si", "diamond", a=5.43)

    # Configure the mock to return a failed process
    mock_subprocess_run.return_value = subprocess.CompletedProcess(
        args=["pw.x", "-in", "qe.in"],
        returncode=1,
        stdout="Starting QE...\nError!",
        stderr="MPI error",
    )

    # Act
    db_id = labelling_engine.execute(atoms)

    # Assert
    retrieved_atoms, result = mock_db.get(db_id)
    assert result.was_successful is False
    assert "QE exited with error code 1" in result.error_message
    assert "MPI error" in result.error_message

@patch("subprocess.run")
def test_labelling_engine_execute_convergence_failure(mock_subprocess_run, labelling_engine: LabellingEngine, mock_db: AseDB):
    """
    Tests the LabellingEngine's handling of a QE run that fails to converge.
    """
    # Arrange
    atoms = bulk("Si", "diamond", a=5.43)
    failure_output = (SAMPLE_DATA_DIR / "failure.out").read_text()

    # Configure the mock for a successful exit but failed convergence in the output
    mock_subprocess_run.return_value = subprocess.CompletedProcess(
        args=["pw.x", "-in", "qe.in"],
        returncode=0,
        stdout=failure_output,
        stderr="",
    )

    # Act
    db_id = labelling_engine.execute(atoms)

    # Assert
    retrieved_atoms, result = mock_db.get(db_id)
    assert result.was_successful is False
    assert result.error_message == "SCF calculation did not converge."
