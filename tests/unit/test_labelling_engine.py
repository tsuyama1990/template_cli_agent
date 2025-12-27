import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path
import numpy as np

from ase.build import bulk

from src.mlip_autopipec.data.database import AseDB
from src.mlip_autopipec.data.models import DFTResult
from src.mlip_autopipec.modules.c_labelling_engine import LabellingEngine

# Sample Quantum Espresso output for a successful run
SUCCESSFUL_QE_OUTPUT = """
!    total energy              =     -11.39867073 Ry
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  Si   force =    -0.00000000    0.00000000    0.00000000
     atom    2 type  Si   force =     0.00000000   -0.00000000    0.00000000

     total   stress  (Ry/bohr**3)                       pbc
      -0.00000012   0.00000000  -0.00000000          ...
       0.00000000  -0.00000012   0.00000000          ...
      -0.00000000   0.00000000  -0.00000012          ...

End of self-consistent calculation

JOB DONE.
"""

# Sample output for a failed (non-converged) run
FAILED_QE_OUTPUT = """
     Self-consistent calculation failed
     c_bands:      2 eigenvalues not converged
     c_bands:      2 eigenvalues not converged
     c_bands:      1 eigenvalues not converged
     c_bands:      1 eigenvalues not converged

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine c_bands (2):
     SCF NOT CONVERGED
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     stopping ...
"""

@pytest.fixture
def mock_db(tmp_path):
    """Fixture to create a mock AseDB instance for testing."""
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)
    # Mock the actual write method to avoid file I/O
    db.write = MagicMock(return_value=1)
    return db

@pytest.fixture
def silicon_atoms():
    """Fixture to create a standard silicon Atoms object."""
    return bulk('Si', 'diamond', a=5.43)

@patch('shutil.which', return_value='/usr/bin/pw.x')
def test_labelling_engine_init_success(mock_which):
    """Test that LabellingEngine initialises correctly if pw.x is found."""
    engine = LabellingEngine(qe_command="pw.x", db=MagicMock())
    assert engine._qe_command == "pw.x"

@patch('shutil.which', return_value=None)
def test_labelling_engine_init_fail(mock_which):
    """Test that LabellingEngine raises FileNotFoundError if pw.x is not found."""
    with pytest.raises(FileNotFoundError):
        LabellingEngine(qe_command="pw.x", db=MagicMock())

@patch('subprocess.run')
@patch('shutil.which', return_value='/usr/bin/pw.x')
def test_execute_successful_run(mock_which, mock_subprocess_run, mock_db, silicon_atoms, tmp_path):
    """Test the full execute workflow for a successful QE run."""
    # Configure the mock to return a successful process result
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = SUCCESSFUL_QE_OUTPUT
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(qe_command="mpirun -np 2 pw.x", db=mock_db)

    db_id = engine.execute(silicon_atoms, work_dir=tmp_path / "qe_temp")

    # Verify that the subprocess was called correctly
    expected_command = ["mpirun", "-np", "2", "pw.x", "-in", str(tmp_path / "qe_temp/qe_input.in")]
    mock_subprocess_run.assert_called_once_with(
        expected_command,
        capture_output=True,
        text=True,
        cwd=tmp_path / "qe_temp"
    )

    # Verify that the database write method was called with the correct data
    mock_db.write.assert_called_once()
    call_args = mock_db.write.call_args
    # The first argument is the atoms object
    assert call_args[0][0] == silicon_atoms
    # The keyword argument is the result object
    result = call_args[1]['result']

    assert result.was_successful is True
    assert result.error_message is None
    assert result.total_energy_ev == pytest.approx(-11.39867073 * 13.6057, abs=1e-6)
    assert np.allclose(result.forces, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

@patch('subprocess.run')
@patch('shutil.which', return_value='/usr/bin/pw.x')
def test_execute_failed_run(mock_which, mock_subprocess_run, mock_db, silicon_atoms, tmp_path):
    """Test the execute workflow for a failed (non-converged) QE run."""
    mock_process = MagicMock()
    mock_process.returncode = 1  # QE often returns a non-zero code on failure
    mock_process.stdout = FAILED_QE_OUTPUT
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(qe_command="pw.x", db=mock_db)

    db_id = engine.execute(silicon_atoms, work_dir=tmp_path / "qe_temp")

    # Verify that the DB was still called, but with failure data
    mock_db.write.assert_called_once()
    result = mock_db.write.call_args[1]['result']

    assert result.was_successful is False
    assert result.total_energy_ev is None
    assert result.forces is None
    assert "SCF did not converge" in result.error_message

@patch('subprocess.run')
@patch('shutil.which', return_value='/usr/bin/pw.x')
def test_execute_subprocess_error(mock_which, mock_subprocess_run, mock_db, silicon_atoms, tmp_path):
    """Test the workflow when subprocess.run itself raises an exception."""
    mock_subprocess_run.side_effect = OSError("Something went wrong")

    engine = LabellingEngine(qe_command="pw.x", db=mock_db)

    db_id = engine.execute(silicon_atoms, work_dir=tmp_path / "qe_temp")

    # Verify DB was called with a generic error message
    mock_db.write.assert_called_once()
    result = mock_db.write.call_args[1]['result']

    assert result.was_successful is False
    assert "Something went wrong" in result.error_message
