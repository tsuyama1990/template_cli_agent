"""Unit tests for the database wrapper."""

from pathlib import Path

import pytest
from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import AseDBWrapper


def test_db_creation_and_write(tmp_path: Path) -> None:
    # Corresponds to SPEC.md, Section 2: Component Blueprint (`storage/database.py`)
    # Verifies that the AseDBWrapper correctly creates a new database file.
    """Test that the database file is created on first write."""
    db_path = tmp_path / "test.db"
    assert not db_path.exists()

    atoms = bulk("Cu", "fcc", a=3.6)
    with AseDBWrapper(db_path) as db_wrapper:
        db_wrapper.write_structures([atoms])

    assert db_path.exists()


def test_write_and_read_integrity(tmp_path: Path) -> None:
    # Corresponds to UAT-C1-004: Verify the Contents and Integrity of the Output Database
    # Verifies that data written to the database can be read back without corruption.
    """Test that writing and reading back an Atoms object preserves its integrity."""
    db_path = tmp_path / "test.db"
    atoms_to_write = bulk("Au", "fcc", a=4.0)

    with AseDBWrapper(db_path) as db_wrapper:
        db_wrapper.write_structures([atoms_to_write])

    db = connect(db_path)
    retrieved_atoms = db.get_atoms(id=1)

    assert len(retrieved_atoms) == len(atoms_to_write)
    assert all(retrieved_atoms.get_atomic_numbers() == atoms_to_write.get_atomic_numbers())
    assert (retrieved_atoms.get_cell() == atoms_to_write.get_cell()).all()
    assert (retrieved_atoms.get_positions() == atoms_to_write.get_positions()).all()


def test_write_outside_context_manager(tmp_path: Path) -> None:
    # Corresponds to SPEC.md, Section 4: Database Wrapper
    # Verifies the robustness of the AseDBWrapper by ensuring it prevents
    # unsafe operations outside of its context manager.
    """Test that writing outside the context manager raises an error."""
    db_path = tmp_path / "test.db"
    db_wrapper = AseDBWrapper(db_path)
    atoms = bulk("Cu", "fcc", a=3.6)
    with pytest.raises(RuntimeError):
        db_wrapper.write_structures([atoms])
