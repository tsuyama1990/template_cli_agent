from pathlib import Path

import pytest
from ase.build import bulk

from mlip_autopipec.domain_models import DFTResult
from mlip_autopipec.utils.ase_db import AseDB


@pytest.fixture
def temp_db_path(tmp_path: Path) -> Path:
    """Creates a temporary database path for testing."""
    return tmp_path / "test.db"


def test_asedb_write_and_get(temp_db_path: Path):
    """Tests writing to and retrieving from the AseDB."""
    asedb = AseDB(temp_db_path)
    atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        energy=-100.0,
        forces=[[0.0, 0.0, 0.0]] * len(atoms),
        stress=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    )

    # Write to the database
    uid = asedb.write(atoms, dft_result, was_successful=True)
    assert isinstance(uid, int)

    # Retrieve from the database
    retrieved_atoms = asedb.get(uid)
    assert retrieved_atoms is not None
    assert len(retrieved_atoms) == len(atoms)
    assert retrieved_atoms.get_potential_energy() == pytest.approx(-100.0)
    assert retrieved_atoms.get_forces().shape == (len(atoms), 3)


def test_asedb_get_non_existent(temp_db_path: Path):
    """Tests that retrieving a non-existent UID returns None."""
    asedb = AseDB(temp_db_path)
    assert asedb.get(999) is None
