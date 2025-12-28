# Description: Unit tests for the AseDB wrapper class.
import os

import pytest
from ase.atoms import Atoms
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


@pytest.fixture
def temp_db_path(tmpdir):
    """Creates a temporary database file for testing."""
    db_path = os.path.join(str(tmpdir), "test.db")
    yield db_path
    if os.path.exists(db_path):
        os.remove(db_path)


def test_asedb_write_and_get_successful_result(temp_db_path):
    """
    Tests writing a successful DFT result to the database and retrieving it.
    """
    db = AseDB(temp_db_path)
    si_atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0, 0.0, 0.0]] * len(si_atoms),
        stress=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        was_successful=True,
    )

    db_id = db.write(si_atoms, dft_result)
    assert isinstance(db_id, int)

    retrieved_atoms, kvp = db.get(db_id)

    assert isinstance(retrieved_atoms, Atoms)
    assert retrieved_atoms.get_potential_energy() == -100.0
    assert kvp["was_successful"] is True
    assert "error_message" not in kvp


def test_asedb_write_and_get_failed_result(temp_db_path):
    """
    Tests writing a failed DFT result to the database and retrieving it.
    """
    db = AseDB(temp_db_path)
    si_atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        total_energy_ev=0.0,
        forces=[],
        stress=[],
        was_successful=False,
        error_message="SCF did not converge",
    )

    db_id = db.write(si_atoms, dft_result)
    retrieved_atoms, kvp = db.get(db_id)

    assert kvp["was_successful"] is False
    assert kvp["error_message"] == "SCF did not converge"
    assert retrieved_atoms.calc is None
