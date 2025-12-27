"""Unit tests for the AseDB interface."""
import pytest
from ase.atoms import Atoms
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
import os
import numpy as np

@pytest.fixture
def db_path(tmp_path):
    """Fixture for a temporary database path."""
    return os.path.join(tmp_path, "test.db")

def test_write_and_get_successful_result(db_path):
    """Test writing and reading a successful DFT result."""
    db = AseDB(db_path)
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
    result = DFTResult(
        total_energy_ev=-10.0,
        forces=[[0, 0, 0.1], [0, 0, -0.1]],
        stress=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        was_successful=True,
    )

    db_id = db.write(atoms, result)
    assert isinstance(db_id, int)

    retrieved_data = db.get(db_id)
    assert retrieved_data is not None
    assert isinstance(retrieved_data["atoms"], Atoms)
    assert retrieved_data["atoms"].get_chemical_symbols() == ['H', 'H']
    assert np.allclose(retrieved_data["total_energy_ev"], -10.0)
    assert np.allclose(retrieved_data["forces"], [[0, 0, 0.1], [0, 0, -0.1]])
    assert np.allclose(retrieved_data["stress"], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert retrieved_data["was_successful"] is True
    assert retrieved_data["error_message"] is None

def test_write_and_get_failed_result(db_path):
    """Test writing and reading a failed DFT result."""
    db = AseDB(db_path)
    atoms = Atoms('O', positions=[[0, 0, 0]])
    result = DFTResult(
        total_energy_ev=0.0,
        forces=[[0, 0, 0]],
        stress=[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        was_successful=False,
        error_message="SCF failed to converge",
    )

    db_id = db.write(atoms, result)
    retrieved_data = db.get(db_id)

    assert retrieved_data is not None
    assert retrieved_data["was_successful"] is False
    assert retrieved_data["error_message"] == "SCF failed to converge"

def test_get_nonexistent_id(db_path):
    """Test getting a non-existent ID returns None."""
    db = AseDB(db_path)
    assert db.get(999) is None
