import pytest
import subprocess
from unittest.mock import patch, MagicMock
from ase.build import bulk
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.database import AseDB

@pytest.fixture
def mock_db(tmp_path):
    db_path = tmp_path / "test.db"
    return AseDB(str(db_path))

@patch('subprocess.run')
def test_labelling_engine_success(mock_subprocess, mock_db):
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock a successful QE run
    mock_process = MagicMock()
    mock_process.stdout = """
         JOB DONE.
         !    total energy              =     -100.0 Ry
         Forces acting on atoms (cartesian axes, Ry/au):

         atom    1 type  1   force =     0.0 0.0 0.0
         atom    2 type  1   force =     0.0 0.0 0.0

         total   stress  (Ry/bohr**3)                       pbar
          -0.1   -0.2   -0.3          -1.4   -2.9   -4.3
          -0.4   -0.5   -0.6          -5.8   -7.2   -8.7
          -0.7   -0.8   -0.9         -10.1  -11.6  -13.0
    """
    mock_process.stderr = ""
    mock_subprocess.return_value = mock_process

    engine = LabellingEngine(
        qe_command="pw.x",
        parameters={"pseudopotentials": {"Si": "Si.upf"}},
        db=mock_db
    )

    db_id = engine.execute(atoms)

    result_metadata = mock_db.read(db_id)
    assert result_metadata['was_successful'] is True

@patch('subprocess.run')
def test_labelling_engine_failure(mock_subprocess, mock_db):
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock a failed QE run
    mock_subprocess.side_effect = subprocess.CalledProcessError(
        1, "pw.x", stderr="Error: SCF not converged"
    )

    engine = LabellingEngine(
        qe_command="pw.x",
        parameters={"pseudopotentials": {"Si": "Si.upf"}},
        db=mock_db
    )

    db_id = engine.execute(atoms)

    result_metadata = mock_db.read(db_id)
    assert result_metadata['was_successful'] is False
    assert "SCF not converged" in result_metadata['error_message']
