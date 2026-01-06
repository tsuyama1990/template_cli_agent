from pathlib import Path

import pytest
from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import AseDBWrapper


def test_asedb_wrapper_write_read(tmp_path: Path) -> None:
    """Test that the AseDBWrapper can write and read structures."""
    db_path = tmp_path / "test.db"
    atoms = bulk("Cu", "fcc", a=3.6)

    # Write to the database
    with AseDBWrapper(db_path) as db:
        db.write_structures([atoms])

    # Read from the database to verify
    with connect(db_path) as db_conn:  # type: ignore[no-untyped-call]
        row = db_conn.get(id=1)
        read_atoms = row.toatoms()

    assert len(read_atoms) == len(atoms)
    assert all(read_atoms.numbers == atoms.numbers)


def test_write_outside_context_manager(tmp_path: Path) -> None:
    """Test that writing outside the context manager raises an error."""
    db_path = tmp_path / "test.db"
    atoms = bulk("Cu", "fcc", a=3.6)
    db_wrapper = AseDBWrapper(db_path)
    with pytest.raises(ConnectionError):
        db_wrapper.write_structures([atoms])
