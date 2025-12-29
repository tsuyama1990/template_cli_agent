import pytest
from unittest.mock import MagicMock, patch, call
from ase import Atoms

from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult

@pytest.fixture
def mock_db():
    """Fixture for a mocked AseDB instance."""
    db = MagicMock(spec=AseDB)
    db.write.return_value = 1  # Assume DB always returns ID 1
    return db

@pytest.fixture
def successful_qe_output():
    """Fixture to provide a sample successful QE output string."""
    with open("tests/unit/test_data/qe_outputs/h_atom_successful_run.out") as f:
        return f.read()

@patch('subprocess.run')
@patch('mlip_autopipec.utils.dft_utils.generate_qe_input')
def test_labelling_engine_execute_successful_run(
    mock_generate_input, mock_subprocess_run, mock_db, successful_qe_output
):
    """
    Tests the LabellingEngine's execute method for a successful QE run.
    It verifies that input is generated, QE is called, output is parsed,
    and the result is written to the database.
    """
    # Arrange
    mock_generate_input.return_value = "dummy qe input"

    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = successful_qe_output
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(qe_command="mpirun -np 4 pw.x", db=mock_db)
    atoms = Atoms('H', positions=[[0, 0, 0]], cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # Act
    db_id = engine.execute(atoms, {}, {}, (1,1,1))

    # Assert
    mock_generate_input.assert_called_once_with(atoms, {}, {}, (1,1,1))
    mock_subprocess_run.assert_called_once_with(
        ["mpirun", "-np", "4", "pw.x"], # Command is now a list
        capture_output=True,
        text=True,
        input="dummy qe input"
    )

    # Check that the DB's write method was called with a valid DFTResult
    assert mock_db.write.call_count == 1
    call_args, _ = mock_db.write.call_args
    assert call_args[0] == atoms
    assert isinstance(call_args[1], DFTResult)
    assert call_args[1].was_successful is True
    assert db_id == 1
