import pytest
from unittest.mock import Mock, patch, call
import subprocess

from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine

# Sample successful QE output
SAMPLE_QE_OUTPUT_SUCCESS = """
!    total energy              =      -100.0 Ry
Forces acting on atoms (cartesian axes, Ry/au):
 atom    1 type  1   force =     0.1      0.2      0.3
total stress  (Ry/bohr**3)                (kbar)     P=      -1.12
  1.0   0.4   0.5
  0.4   2.0   0.6
  0.5   0.6   3.0
"""

@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper."""
    db = Mock(spec=AseDBWrapper)

    # Mock row object that can be converted to an Atoms object
    mock_row = Mock()
    mock_row.id = 1
    mock_row.toatoms.return_value = Mock() # A mock atoms object

    db.get_rows_to_label.return_value = [mock_row]
    return db

@pytest.fixture
def dft_config():
    """Fixture for a sample DFTCompute config."""
    return DFTCompute(
        code="quantum_espresso",
        command="mpirun -np 4 pw.x",
        pseudopotentials="SSSP",
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=3.0,
    )

@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
@patch('mlip_autopipec.modules.c_labeling_engine.parse_qe_output')
def test_labeling_engine_success(
    mock_parse, mock_create, mock_run, mock_open, mock_makedirs,
    dft_config, mock_db_wrapper
):
    """Test the LabelingEngine's successful execution path."""
    # Arrange
    mock_create.return_value = "dummy_qe_input"
    mock_parse.return_value = {'energy': -1360.5, 'forces': [[1,1,1]], 'stress': [1,2,3,4,5,6]}

    # Mock the result of subprocess.run
    mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0)

    # Act
    engine = LabelingEngine(dft_config, mock_db_wrapper)
    engine.execute()

    # Assert
    mock_db_wrapper.get_rows_to_label.assert_called_once()
    mock_create.assert_called_once()

    # Check that subprocess.run was called correctly
    expected_command = ["mpirun", "-np", "4", "pw.x", "-in", "calc_1/qe.in"]
    mock_run.assert_called_once()
    called_command = mock_run.call_args[0][0]
    assert called_command == expected_command

    mock_parse.assert_called_once()
    mock_db_wrapper.update_row_with_dft_results.assert_called_once_with(1, mock_parse.return_value)

@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
def test_labeling_engine_dft_failure(
    mock_create, mock_run, mock_open, mock_makedirs, dft_config, mock_db_wrapper
):
    """Test how the LabelingEngine handles a failed DFT calculation."""
    # Arrange
    mock_run.side_effect = subprocess.CalledProcessError(returncode=1, cmd="pw.x", stderr=b"QE error")

    # Act
    engine = LabelingEngine(dft_config, mock_db_wrapper)
    engine.execute()

    # Assert
    mock_run.assert_called_once()
    mock_db_wrapper.update_row_with_dft_results.assert_not_called()

@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.parse_qe_output')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
def test_labeling_engine_parsing_failure(
    mock_create, mock_parse, mock_run, mock_open, mock_makedirs, dft_config, mock_db_wrapper
):
    """Test how the LabelingEngine handles a failure in parsing the output."""
    # Arrange
    mock_parse.return_value = None # Simulate parsing failure
    mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0)

    # Act
    engine = LabelingEngine(dft_config, mock_db_wrapper)
    engine.execute()

    # Assert
    mock_run.assert_called_once()
    mock_parse.assert_called_once()
    mock_db_wrapper.update_row_with_dft_results.assert_not_called()
