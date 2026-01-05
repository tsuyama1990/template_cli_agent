"""Wrapper for interacting with the ASE database."""
from pathlib import Path
from typing import List, Optional

from ase.atoms import Atoms
from ase.db import connect
from ase.db.core import Database


class AseDBWrapper:
    """A context manager for the ASE database."""

    def __init__(self, db_path: Path):
        self.db_path = db_path
        self._connection: Optional[Database] = None

    def __enter__(self) -> "AseDBWrapper":
        """Open the database connection."""
        self._connection = connect(self.db_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close the database connection."""
        # The ASE DB object handles its own connection closing via its destructor.
        # We just need to clear our reference to indicate the context is closed.
        self._connection = None

    def write_structures(self, structures: List[Atoms]):
        """Write a list of Atoms objects to the database."""
        if self._connection is None:
            raise RuntimeError(
                "Database connection is not open. Use within a 'with' statement."
            )
        for atoms in structures:
            self._connection.write(atoms)
