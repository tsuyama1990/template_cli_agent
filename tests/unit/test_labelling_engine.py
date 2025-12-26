import pytest
from unittest.mock import patch, MagicMock
from ase.build import bulk

from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult

@pytest.fixture
def mock_db(tmp_path):
    """A fixture to create a mock AseDB for testing."""
    db_path = tmp_path / "test.db"
    return AseDB(db_path)

def test_labelling_engine_success(mock_db):
    """
    Tests the LabellingEngine for a successful QE run.
    """
    engine = LabellingEngine("pw.x", mock_db)
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock a successful subprocess run
    mock_run = MagicMock()
    mock_run.returncode = 0
    mock_run.stdout = """
    !    total energy              =     -150.0 Ry
    Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =    -0.00000001   -0.00000002    0.00000003
     atom    2 type  1   force =     0.00000001    0.00000002   -0.00000003
    """

    with patch("subprocess.run", return_value=mock_run):
        db_id = engine.execute(atoms)

    assert isinstance(db_id, int)

    # Verify the data was written correctly to the DB
    retrieved_data = DFTResult(**mock_db.get_row(db_id))
    assert retrieved_data.was_successful is True
    assert retrieved_data.total_energy_ev == pytest.approx(-150.0 * 13.6057)

def test_labelling_engine_qe_failure(mock_db):
    """
    Tests the LabellingEngine for a failed QE run (non-zero exit code).
    """
    engine = LabellingEngine("pw.x", mock_db)
    atoms = bulk("Si", "diamond", a=5.43)

    # Mock a failed subprocess run
    mock_run = MagicMock()
    mock_run.returncode = 1
    mock_run.stderr = "Quantum Espresso crashed."

    with patch("subprocess.run", return_value=mock_run):
        db_id = engine.execute(atoms)

    retrieved_data = DFTResult(**mock_db.get_row(db_id))
    assert retrieved_data.was_successful is False
    assert "Quantum Espresso crashed" in retrieved_data.error_message
