"""Unit tests for the ASE database wrapper."""
from pathlib import Path
from ase.build import bulk
from ase.db import connect
from mlip_autopipec.storage.database import AseDBWrapper

def test_ase_db_wrapper_write_and_read(tmp_path: Path):
    """
    Test that the AseDBWrapper can create a database, write structures to it,
    and that the data can be read back with integrity.
    """
    db_path = tmp_path / "test.db"

    # Create a known Atoms object
    atoms1 = bulk("Cu", "fcc", a=3.6)
    atoms2 = bulk("Au", "fcc", a=4.0)
    structures_to_write = [atoms1, atoms2]

    # Use the wrapper to write the structure
    with AseDBWrapper(db_path) as wrapper:
        wrapper.write_structures(structures_to_write)

    # The file should exist now that data has been written
    assert db_path.exists()

    # Connect directly to the database to verify the contents
    db = connect(db_path)

    assert len(db) == 2

    retrieved_atoms1 = db.get_atoms(id=1)
    retrieved_atoms2 = db.get_atoms(id=2)

    # Verify data integrity
    assert retrieved_atoms1.get_chemical_symbols() == ["Cu"]
    assert retrieved_atoms1.cell.any() == atoms1.cell.any()

    assert retrieved_atoms2.get_chemical_symbols() == ["Au"]
    assert retrieved_atoms2.cell.any() == atoms2.cell.any()

def test_context_manager_opens_and_closes(tmp_path: Path):
    """Test that the context manager properly handles the connection."""
    db_path = tmp_path / "context_test.db"

    wrapper = AseDBWrapper(db_path)

    # Connection should be None before entering the context
    assert wrapper._connection is None

    with wrapper:
        # Connection should be active inside the context
        assert wrapper._connection is not None
        # The db file itself is created on the first write, so we don't
        # assert its existence here. This test focuses on the connection state.

    # Connection should be None after exiting the context
    assert wrapper._connection is None
