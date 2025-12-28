import subprocess
from unittest.mock import MagicMock, mock_open, patch

import pytest
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

# --- Test Data ---
MOCK_SUCCESS_OUTPUT = """
!    total energy              =      -15.83228912 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000000    -0.00000000    -0.00000000
     atom    2 type  1   force =     0.00000000     0.00000000     0.00000000

total stress  (Ry/bohr**3)  (kbar)     P=      -2.67
  -0.00016867   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000
"""

MOCK_FAILURE_OUTPUT = "SCF NOT CONVERGED"


@pytest.fixture
def silicon_atoms():
    """Provides a silicon Atoms object for testing."""
    return bulk("Si", "diamond", a=5.43)


def test_generate_qe_input(silicon_atoms):
    """Tests that the QE input generation is syntactically plausible."""
    params = {"ecutwfc": 70, "k_points": [4, 4, 4], "pseudo_dir": "test_pseudos"}
    input_str = generate_qe_input(silicon_atoms, params)
    assert "&CONTROL" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 70" in input_str
    assert "ATOMIC_POSITIONS" in input_str
    # Corrected assertion to match the actual output of ase.build.bulk
    assert "  Si 0.00000000 0.00000000 0.00000000" in input_str
    assert "  Si 1.35750000 1.35750000 1.35750000" in input_str
    assert "CELL_PARAMETERS" in input_str
    assert "4 4 4 0 0 0" in input_str
    assert "pseudo_dir = 'test_pseudos'" in input_str


def test_parse_qe_output_success():
    """Tests parsing a successful QE output."""
    energy, forces, stress, success, error = parse_qe_output(MOCK_SUCCESS_OUTPUT)
    assert success is True
    assert error is None
    assert energy == pytest.approx(-15.83228912 * 13.605693, abs=1e-4)
    assert len(forces) == 2
    assert forces[0][0] == pytest.approx(0.0)
    assert stress is not None
    assert stress[0][0] == pytest.approx(-0.00016867 * (1 / 160.21766), abs=1e-6)


def test_parse_qe_output_failure():
    """Tests parsing a failed QE output."""
    energy, forces, stress, success, error = parse_qe_output(MOCK_FAILURE_OUTPUT)
    assert success is False
    assert "Could not find total energy" in error


@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
@patch("builtins.open", new_callable=mock_open, read_data=MOCK_SUCCESS_OUTPUT)
def test_labelling_engine_execute_success(mock_file, mock_subprocess_run, silicon_atoms):
    """Tests the successful execution workflow of the LabellingEngine."""
    mock_db = MagicMock(spec=AseDB)
    mock_db.write.return_value = 1
    mock_subprocess_run.return_value = MagicMock(returncode=0)

    engine = LabellingEngine(qe_command="pw.x", db=mock_db, default_params={})
    db_id = engine.execute(silicon_atoms)

    assert db_id == 1
    mock_subprocess_run.assert_called_once()
    mock_db.write.assert_called_once()

    # Check the arguments passed to the mock db.write method
    args, _ = mock_db.write.call_args
    # args[0] is the atoms object, args[1] is the dictionary of results
    result_dict = args[1]
    result = DFTResult(**result_dict)
    assert result.was_successful is True
    assert result.total_energy_ev < 0


@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
def test_labelling_engine_execute_subprocess_error(mock_subprocess_run, silicon_atoms):
    """Tests the LabellingEngine's handling of a subprocess failure."""
    mock_db = MagicMock(spec=AseDB)
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, "cmd", stderr="SCF failed")

    engine = LabellingEngine(qe_command="pw.x", db=mock_db, default_params={})
    engine.execute(silicon_atoms)

    mock_db.write.assert_called_once()
    args, _ = mock_db.write.call_args
    result_dict = args[1]
    result = DFTResult(**result_dict)
    assert result.was_successful is False
    assert "execution failed" in result.error_message
