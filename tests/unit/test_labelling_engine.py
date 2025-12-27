import pytest
from unittest.mock import patch, MagicMock
from ase.build import bulk
from src.mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from src.mlip_autopipec.data.database import AseDB

# Re-using the sample output from the dft_utils test
SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =    -111.92138384 Ry
     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =    -0.00000002    -0.00000002    -0.00000002
     atom    2 type  1   force =     0.00000002     0.00000002     0.00000002

     Total Force =     0.00000004     Total SCF correction =     0.00000000

          total stress  (Ry/bohr**3)            (kbar)       P=   -2.29
   -0.00001452   0.00000000   0.00000000      -0.23      0.00      0.00
    0.00000000  -0.00001452   0.00000000       0.00     -0.23      0.00
    0.00000000   0.00000000  -0.00001452       0.00      0.00     -0.23
"""

SAMPLE_FAILURE_OUTPUT = """
     Maximum number of SCF cycles reached.
     Message from routine cdiaghg:
     not converged in 100 iterations
"""

@pytest.fixture
def mock_db(tmp_path):
    """Fixture for a mocked AseDB."""
    db_path = tmp_path / "test.db"
    db = AseDB(str(db_path))
    db.write = MagicMock(return_value=1) # Mock the write method
    return db

@patch("subprocess.run")
def test_labelling_engine_success(mock_subprocess_run, mock_db):
    """
    Tests the LabellingEngine for a successful QE run.
    """
    # Mock the subprocess call to return a successful run with sample output
    mock_subprocess_run.return_value = MagicMock(
        returncode=0, stdout=SAMPLE_SUCCESS_OUTPUT, stderr=""
    )

    engine = LabellingEngine(
        qe_command="mpirun -np 4 pw.x",
        db=mock_db,
        pseudo_dir="/pseudos",
        pseudopotentials={"Si": "Si.upf"},
    )

    atoms = bulk("Si", "diamond", a=5.43)
    db_id = engine.execute(atoms)

    # Verify that subprocess.run was called correctly
    assert mock_subprocess_run.called
    args, kwargs = mock_subprocess_run.call_args
    command_list = args[0]
    assert command_list[0:4] == ["mpirun", "-np", "4", "pw.x"]

    # Verify that the database write method was called
    assert mock_db.write.called
    call_args, call_kwargs = mock_db.write.call_args
    written_atoms = call_args[0]
    written_result = call_kwargs['result']

    assert written_atoms == atoms
    assert written_result.was_successful
    assert written_result.total_energy_ev == pytest.approx(-111.92138384 * 13.6057)
    assert db_id == 1


@patch("subprocess.run")
def test_labelling_engine_failure(mock_subprocess_run, mock_db):
    """
    Tests the LabellingEngine for a failed QE run.
    """
    mock_subprocess_run.return_value = MagicMock(
        returncode=1, stdout=SAMPLE_FAILURE_OUTPUT, stderr="Error"
    )

    engine = LabellingEngine(
        qe_command="pw.x",
        db=mock_db,
        pseudo_dir="/pseudos",
        pseudopotentials={"Si": "Si.upf"},
    )

    atoms = bulk("Si", "diamond", a=5.43)
    db_id = engine.execute(atoms)

    # Verify that the DB write method was called with a failure status
    assert mock_db.write.called
    _, write_kwargs = mock_db.write.call_args
    written_result = write_kwargs['result']

    assert not written_result.was_successful
    assert "SCF did not converge" in written_result.error_message
    assert db_id == 1
