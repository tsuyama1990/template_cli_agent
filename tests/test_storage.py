from pathlib import Path

import pytest
from ase.build import bulk
from ase.db import connect

from mlip_autopipec.storage.database import AseDBWrapper


def test_ase_db_wrapper(tmp_path: Path):
    """Tests that the AseDBWrapper can write and read structures."""
    db_path = tmp_path / "test.db"
    wrapper = AseDBWrapper(db_path)

    # Create a test structure
    atoms = bulk("Cu", "fcc", a=3.6)

    # Write the structure to the database
    wrapper.write_structures([atoms])

    # Read the structure back from the database
    with connect(db_path) as db:
        read_atoms = db.get_atoms(id=1)

    # Check that the read structure is the same as the original
    assert len(read_atoms) == len(atoms)
    assert (read_atoms.get_atomic_numbers() == atoms.get_atomic_numbers()).all()
    assert (read_atoms.get_positions() == atoms.get_positions()).all()
