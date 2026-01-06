# tests/mlip_autopipec/unit/test_storage.py
"""Unit tests for the database manager."""
import tempfile
import os
from ase import Atoms
from ase.db import connect
from ase.calculators.singlepoint import SinglePointCalculator
from mlip_autopipec.storage.database_manager import DatabaseManager


def test_database_manager_writes_structures() -> None:
    """Test that the DatabaseManager correctly writes a list of Atoms to a file."""
    # 1. Create a temporary file for the database
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        db_path = tmp.name

    # 2. Create mock Atoms objects with calculators
    atoms1 = Atoms("H")
    atoms1.calc = SinglePointCalculator(atoms1, energy=-1.0)
    atoms2 = Atoms("He")
    atoms2.calc = SinglePointCalculator(atoms2, energy=-2.0)
    structures = [atoms1, atoms2]

    # 3. Initialize the manager and write the structures
    manager = DatabaseManager(db_path)
    manager.write_structures(structures)

    # 4. Assertions
    # Connect to the database and verify the contents
    with connect(db_path) as db:
        assert db.count() == 2
        row1 = db.get(id=1)
        row2 = db.get(id=2)
        assert row1.symbols == ["H"]
        assert row2.energy == -2.0

    # Clean up the temporary file
    os.remove(db_path)


def test_database_manager_handles_empty_list() -> None:
    """Test that the DatabaseManager handles an empty list of structures gracefully."""
    # 1. Create a temporary file
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        db_path = tmp.name

    # 2. Initialize the manager and write an empty list
    manager = DatabaseManager(db_path)
    manager.write_structures([])

    # 3. Assertions
    with connect(db_path) as db:
        assert db.count() == 0

    # Clean up
    os.remove(db_path)
