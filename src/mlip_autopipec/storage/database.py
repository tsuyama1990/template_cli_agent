"""ASE database wrapper."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from ase.db import connect

if TYPE_CHECKING:
    from ase import Atoms
    from ase.db.core import Connection


class AseDBWrapper:
    """A context manager for writing to an ASE database."""

    def __init__(self, db_path: str | Path) -> None:
        """
        Initialise the AseDBWrapper.
        Args:
            db_path: The path to the database file.
        """
        self.db_path = db_path
        self._connection: Connection | None = None

    def __enter__(self) -> AseDBWrapper:
        """
        Enter the context manager and connect to the database.
        """
        self._connection = connect(self.db_path)
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """
        Exit the context manager.
        The ASE db connection is managed implicitly and doesn't need explicit closing.
        """
        self._connection = None

    def write_structures(self, structures: list[Atoms]) -> None:
        """
        Write a list of Atoms objects to the database.
        Args:
            structures: A list of ase.Atoms objects to write.
        """
        if self._connection is None:
            raise RuntimeError("Database connection not open. Use within a 'with' statement.")

        for atoms in structures:
            self._connection.write(atoms)
