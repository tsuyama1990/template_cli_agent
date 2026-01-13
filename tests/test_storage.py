"""Unit tests for the storage component."""
from pathlib import Path

import pytest
from ase.build import bulk

from mlip_autopipec.storage.database import AseDBWrapper


def test_db_wrapper_creation_and_write(tmp_path: Path) -> None:
    """Test database creation and writing structures."""
    db_path = tmp_path / "test.db"
    atoms1 = bulk("Cu", "fcc", a=3.6)
    atoms2 = bulk("Au", "fcc", a=4.0)
    structures = [atoms1, atoms2]

    assert not db_path.exists()

    with AseDBWrapper(db_path) as db:
        db.write_structures(structures)

    assert db_path.exists()

    # Verify content with ase.db
    from ase.db import connect

    db_conn = connect(db_path)
    assert len(db_conn) == 2
    retrieved_atoms1 = db_conn.get_atoms(id=1)
    retrieved_atoms2 = db_conn.get_atoms(id=2)
    assert retrieved_atoms1.get_chemical_symbols() == ["Cu"]
    assert retrieved_atoms2.get_chemical_symbols() == ["Au"]


def test_db_wrapper_append(tmp_path: Path) -> None:
    """Test that writing to an existing database appends structures."""
    db_path = tmp_path / "test_append.db"
    atoms1 = bulk("Cu", "fcc", a=3.6)
    atoms2 = bulk("Au", "fcc", a=4.0)

    with AseDBWrapper(db_path) as db:
        db.write_structures([atoms1])

    with AseDBWrapper(db_path) as db:
        db.write_structures([atoms2])

    from ase.db import connect

    db_conn = connect(db_path)
    assert len(db_conn) == 2
