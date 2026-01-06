# -*- coding: utf-8 -*-
"""Unit tests for the AseDBWrapper."""
from pathlib import Path

import pytest
from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import AseDBWrapper


def test_db_wrapper_creates_file_and_writes(tmp_path: Path) -> None:
    """Test that the AseDBWrapper creates a database file and writes to it correctly."""
    db_path = tmp_path / "test.db"
    atoms = bulk("Cu", "fcc", a=3.6)

    assert not db_path.exists()

    with AseDBWrapper(db_path) as db_wrapper:
        db_wrapper.write_structures([atoms])

    # ASE DB is created lazily, so it only exists after a write
    assert db_path.exists()

    # Verify content by reading it back with ASE's native connect
    db = connect(db_path)  # type: ignore
    retrieved_atoms = db.get_atoms(id=1)
    assert len(db) == 1
    assert retrieved_atoms.get_chemical_symbols() == ["Cu"]


def test_db_wrapper_context_manager_behavior(tmp_path: Path) -> None:
    """Test that the wrapper cannot be used outside the 'with' context."""
    db_path = tmp_path / "test.db"
    atoms = bulk("Au", "fcc", a=4.0)
    db_wrapper = AseDBWrapper(db_path)

    # Should raise an error because we are not in the context
    with pytest.raises(RuntimeError, match="Database not connected"):
        db_wrapper.write_structures([atoms])

    # Should work correctly inside the context
    with db_wrapper:
        db_wrapper.write_structures([atoms])

    # Should fail again after exiting the context
    with pytest.raises(RuntimeError, match="Database not connected"):
        db_wrapper.write_structures([atoms])

    # Verify that the write within the context was successful
    db = connect(db_path)  # type: ignore
    assert len(db) == 1
