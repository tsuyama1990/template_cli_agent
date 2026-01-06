# tests/unit/storage/test_database_manager.py
import os
import tempfile
from pathlib import Path
from typing import Generator

import pytest
from ase import Atoms
from ase.db import connect

from mlip_autopipec.storage.database_manager import DatabaseManager


@pytest.fixture
def temp_db_path() -> Generator[Path, None, None]:
    """Creates a temporary database file for a test and cleans up after."""
    # Using NamedTemporaryFile to ensure a unique filename
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        path = Path(tmp.name)
    yield path
    # Cleanup: remove the file after the test is done
    if path.exists():
        path.unlink()


def test_database_manager_writes_structures(temp_db_path: Path) -> None:
    """Test that the DatabaseManager can write a list of Atoms objects."""
    db_manager = DatabaseManager()
    structures = [
        Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]]),
        Atoms("O2", positions=[[0, 0, 0], [0, 0, 1.2]]),
    ]

    # Connect to the temporary database and write structures
    with db_manager.connect(str(temp_db_path)) as connection:
        db_manager.write_structures(connection, structures)

    # Verify the contents of the database
    # Re-connect directly to read back the data and assert
    read_con = connect(temp_db_path)  # type: ignore[no-untyped-call]
    assert len(read_con) == 2

    row1 = read_con.get(id=1)
    assert row1.symbols[0] == "H"
    assert row1.natoms == 2

    row2 = read_con.get(id=2)
    assert row2.symbols[0] == "O"
