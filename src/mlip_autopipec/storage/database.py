"""Wrapper for interacting with the ASE database."""
from pathlib import Path
from typing import Sequence

from ase.atoms import Atoms
from ase.db import connect


class AseDBWrapper:
    """A context manager for handling ASE database connections."""

    def __init__(self, db_path: Path | str):
        """
        Initialize the wrapper with the path to the database.

        Args:
            db_path: Path to the ASE database file.
        """
        self.db_path = db_path
        self._db = None

    def __enter__(self):
        """Open the database connection."""
        self._db = connect(self.db_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close the database connection."""
        # The ase.db.connect object doesn't need explicit closing.
        pass

    def write_structures(self, structures: Sequence[Atoms]) -> None:
        """
        Write a sequence of Atoms objects to the database.

        Args:
            structures: A list or sequence of ase.Atoms objects.
        """
        if self._db is None:
            raise ConnectionError("Database is not connected. Use 'with' statement.")

        for atoms in structures:
            self._db.write(atoms)
