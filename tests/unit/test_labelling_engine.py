from unittest.mock import MagicMock, patch

import pytest
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

# Sample successful QE output
SAMPLE_QE_SUCCESS_OUTPUT = """
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000000    -0.00000000    -0.00000000
     atom    2 type  1   force =     0.00000000     0.00000000     0.00000000

     Total force =     0.000000     Total SCF correction =     0.000000

     total stress  (Ry/bohr**3)                (kbar)     P=  -13.71
     -0.00086555   0.00000000   0.00000000        -12.72    0.00    0.00
      0.00000000  -0.00086555   0.00000000          0.00  -12.72    0.00
      0.00000000   0.00000000  -0.00086555          0.00    0.00  -12.72

!    total energy              =     -11.43845671 Ry
     convergence has been achieved in  8 iterations

     JOB DONE.
"""

# Sample failed (non-converged) QE output
SAMPLE_QE_FAIL_OUTPUT = """
     iteration #  8     ecut=    40.00 Ry     beta= 0.70
     total energy =     -11.43845671 Ry
     negative rho (up, down):  2.3E-05 2.3E-05

     convergence NOT achieved after 100 iterations: stopping

     JOB DONE.
"""


@pytest.fixture
def si_atoms():
    """Returns a silicon bulk Atoms object."""
    return bulk("Si", "diamond", a=5.43)


def test_generate_qe_input(si_atoms):
    """Test the generation of a QE input file."""
    params = {
        "control": {"pseudo_dir": "/path/to/pseudos"},
        "system": {"ecutwfc": 60, "ecutrho": 240},
        "electrons": {"conv_thr": 1.0e-8},
        "k_points": [4, 4, 4, 0, 0, 0],
    }
    pseudos = {"Si": "Si.upf"}

    content = generate_qe_input(si_atoms, params, pseudos)

    assert "&CONTROL" in content
    assert "pseudo_dir = '/path/to/pseudos'" in content
    assert "&SYSTEM" in content
    assert "ecutwfc = 60" in content
    assert "nat = 2" in content  # 2 atoms in conventional Si cell
    assert "ATOMIC_SPECIES" in content
    assert "Si " in content
    assert "ATOMIC_POSITIONS {angstrom}" in content
    assert "CELL_PARAMETERS {angstrom}" in content


def test_parse_qe_output_success():
    """Test parsing a successful QE output."""
    success, msg, results = parse_qe_output(SAMPLE_QE_SUCCESS_OUTPUT)

    assert success is True
    assert msg is None
    assert "total_energy_ev" in results
    assert "forces" in results
    assert "stress" in results
    assert abs(results["total_energy_ev"] - (-11.43845671 * 13.605693)) < 1e-4
    assert len(results["forces"]) == 2
    assert len(results["stress"]) == 3
    assert abs(results["forces"][0][0]) < 1e-8


def test_parse_qe_output_failure():
    """Test parsing a failed QE output."""
    success, msg, results = parse_qe_output(SAMPLE_QE_FAIL_OUTPUT)

    assert success is False
    assert "SCF calculation did not converge" in msg
    assert results is None


@patch("src.mlip_autopipec.modules.c_labelling_engine.subprocess.run")
def test_labelling_engine_execute(mock_subprocess_run, si_atoms, tmp_path):
    """Test the LabellingEngine's execute method."""
    mock_db = MagicMock(spec=AseDB)
    mock_db.write.return_value = 1

    # Mock the output file content
    output_file_mock = tmp_path / "qe.out"
    output_file_mock.write_text(SAMPLE_QE_SUCCESS_OUTPUT)

    # Need to patch the tempfile creation to control the output path
    with patch(
        "src.mlip_autopipec.modules.c_labelling_engine.tempfile.TemporaryDirectory"
    ) as mock_tmpdir:
        mock_tmpdir.return_value.__enter__.return_value = str(tmp_path)

        engine = LabellingEngine(
            qe_command="pw.x", db=mock_db, parameters={}, pseudopotentials={"Si": "Si.upf"}
        )

        db_id = engine.execute(si_atoms)

        assert db_id == 1
        mock_subprocess_run.assert_called_once()
        # Check that the command was constructed correctly
        call_args = mock_subprocess_run.call_args[0][0]
        assert "pw.x -in" in call_args
        assert f" > {output_file_mock}" in call_args

        # Check that db.write was called with a successful result
        write_call_args = mock_db.write.call_args[0]
        assert write_call_args[1].was_successful is True
        assert abs(write_call_args[1].total_energy_ev - (-11.43845671 * 13.605693)) < 1e-4
