"""Unit tests for the AseDBWrapper."""
from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import AseDBWrapper


def test_database_creation_on_write(tmp_path):
    """Test that the database is created upon first write."""
    db_path = tmp_path / "test.db"
    wrapper = AseDBWrapper(db_path)
    atoms = bulk("Cu", "fcc", a=3.6)

    assert not db_path.exists()
    wrapper.write_structures([atoms])
    assert db_path.exists()


def test_write_and_read_integrity(tmp_path):
    """Test that writing and reading structures maintains data integrity."""
    db_path = tmp_path / "test.db"
    wrapper = AseDBWrapper(db_path)

    atoms1 = bulk("Cu", "fcc", a=3.6)
    atoms2 = bulk("Au", "fcc", a=4.0)
    structures = [atoms1, atoms2]

    wrapper.write_structures(structures)

    db = connect(db_path)
    assert len(db) == 2

    retrieved_atoms1 = db.get_atoms(id=1)
    retrieved_atoms2 = db.get_atoms(id=2)

    assert retrieved_atoms1.get_chemical_symbols() == ["Cu"]
    assert retrieved_atoms2.get_chemical_symbols() == ["Au"]
