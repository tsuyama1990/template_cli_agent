"""Unit tests for the ASE database wrapper."""

from pathlib import Path

from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import write_structures


def test_ase_db_wrapper(tmp_path: Path):
    """
    Test that the ase_db_wrapper can create a database and write structures to it.
    """
    db_path = tmp_path / "test.db"
    atoms = bulk("Cu", "fcc", a=3.6)
    structures = [atoms]

    write_structures(db_path, structures)

    assert db_path.exists()

    # Verify the contents of the database
    with connect(db_path) as db:
        assert len(db) == 1
        retrieved_atoms = db.get_atoms(id=1)
        assert retrieved_atoms.get_chemical_symbols() == ["Cu"]
