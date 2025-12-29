from unittest.mock import Mock, patch

import pytest
from ase.build import bulk

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine


# Mock AseDB to avoid actual database interaction
@pytest.fixture
def mock_db():
    db = Mock()
    db.write.return_value = 1  # Always return DB ID 1
    return db

# Reusable fixtures from the utils test
@pytest.fixture
def si_atoms():
    return bulk("Si", "diamond", a=5.43)

@pytest.fixture
def qe_parameters():
    return {
        "pseudopotentials": {"Si": "Si.UPF"},
        "k_points": [4, 4, 4],
        "ecutwfc": 60,
        "ecutrho": 240,
    }

# Successful DFT result for mocking parse_qe_output
@pytest.fixture
def successful_dft_result():
    return DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0, 0.0, 0.0]]*2,
        stress=[[0.1]*3]*3,
        was_successful=True
    )

def test_labelling_engine_initialization(mock_db, qe_parameters):
    """
    Tests that the LabellingEngine initializes correctly.
    """
    engine = LabellingEngine(qe_command="pw.x", parameters=qe_parameters, db=mock_db)
    assert engine._qe_command == ["pw.x"]

    # Test that an empty command raises an error
    with pytest.raises(ValueError, match="Quantum Espresso command cannot be empty."):
        LabellingEngine(qe_command="", parameters=qe_parameters, db=mock_db)


@patch("mlip_autopipec.modules.c_labelling_engine.dft_utils")
@patch("subprocess.run")
def test_labelling_engine_execute_successful_run(
    mock_subprocess_run, mock_dft_utils, si_atoms, qe_parameters, mock_db, successful_dft_result
):
    """
    Tests the full execution flow of the LabellingEngine for a successful DFT run.
    Mocks are used to simulate file I/O, subprocess execution, and DB writing.
    """
    # 1. Setup Mocks
    # Mock the return value of the QE process
    mock_process = Mock()
    mock_process.stdout = "Successful QE output"
    mock_process.stderr = ""
    mock_subprocess_run.return_value = mock_process

    # Mock the return value of our parser
    mock_dft_utils.parse_qe_output.return_value = successful_dft_result

    # Mock the input file generator to avoid file system dependency in this test
    mock_dft_utils.generate_qe_input.return_value = "Fake QE input"

    # 2. Execute
    engine = LabellingEngine(qe_command="mpirun -np 2 pw.x", parameters=qe_parameters, db=mock_db)
    db_id = engine.execute(si_atoms)

    # 3. Verify
    # Check that the QE command was called correctly
    mock_subprocess_run.assert_called_once()
    args, kwargs = mock_subprocess_run.call_args
    assert args[0][:4] == ["mpirun", "-np", "2", "pw.x"]
    assert "-in" in args[0]

    # Check that the parser was called with the stdout of the process
    mock_dft_utils.parse_qe_output.assert_called_once_with("Successful QE output")

    # Check that the result was written to the database
    mock_db.write.assert_called_once()
    called_atoms, called_result = mock_db.write.call_args[0]
    assert called_atoms is si_atoms
    assert called_result is successful_dft_result

    # Check that the returned DB ID is correct
    assert db_id == 1

@patch("mlip_autopipec.modules.c_labelling_engine.dft_utils")
@patch("subprocess.run")
def test_labelling_engine_handles_process_error(
    mock_subprocess_run, mock_dft_utils, si_atoms, qe_parameters, mock_db
):
    """
    Tests that the engine correctly handles a failed QE process (e.g., non-zero exit code).
    """
    # 1. Setup Mocks
    mock_process = Mock()
    mock_process.stdout = "Error output"
    mock_process.stderr = "Something went wrong"
    mock_subprocess_run.return_value = mock_process

    # Simulate the parser finding a job failure
    failed_result = DFTResult(was_successful=False, error_message="Job did not finish", total_energy_ev=0, forces=[], stress=[])
    mock_dft_utils.parse_qe_output.return_value = failed_result

    # FIX: Ensure the mocked generator returns a string
    mock_dft_utils.generate_qe_input.return_value = "Fake QE input for failure test"

    # 2. Execute
    engine = LabellingEngine(qe_command="pw.x", parameters=qe_parameters, db=mock_db)
    db_id = engine.execute(si_atoms)

    # 3. Verify
    # Ensure that even with a process error, the parser is called and the result is written
    mock_dft_utils.parse_qe_output.assert_called_once_with("Error output")
    mock_db.write.assert_called_once_with(si_atoms, failed_result)
    assert db_id == 1
