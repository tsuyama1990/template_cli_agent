"""Provides a wrapper for interacting with an ASE (Atomic Simulation Environment) database.

This module encapsulates all database-related operations, offering a clean,
high-level API to the rest of the application. This decouples the core logic
from the specifics of the storage backend, making the system more modular and
easier to maintain. The primary component is the `AseDBWrapper`, which acts as a
context manager for safe and reliable database connections.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from ase.db import connect

if TYPE_CHECKING:
    from ase import Atoms
    from ase.db.core import Connection


class AseDBWrapper:
    """A context manager for safely writing to an ASE database.

    This class handles the connection and management of the output SQLite
    database file used by ASE. By implementing the context manager protocol,
    it ensures that database connections are handled correctly, even in the
    event of errors.

    Attributes:
        db_path (str | Path): The file path for the database.

    """

    def __init__(self, db_path: str | Path) -> None:
        """Initialise the AseDBWrapper.

        Args:
            db_path: The path to the database file.

        """
        self.db_path = db_path
        self._connection: Connection | None = None

    def __enter__(self) -> AseDBWrapper:
        """Enter the context manager and connect to the database.

        Returns:
            The instance of the AseDBWrapper itself.

        """
        self._connection = connect(self.db_path)
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Exit the context manager.

        The ASE db connection is managed implicitly and doesn't need explicit
        closing, but we reset the internal connection state for safety.

        Args:
            exc_type: The exception type if an exception was raised.
            exc_val: The exception value.
            exc_tb: The traceback.

        """
        self._connection = None

    def write_structures(self, structures: list[Atoms]) -> None:
        """Write a list of Atoms objects to the database.

        Args:
            structures: A list of `ase.Atoms` objects to write.

        Raises:
            RuntimeError: If the method is called outside of a `with` block.

        """
        if self._connection is None:
            raise RuntimeError("Database connection not open. Use within a 'with' statement.")

        for atoms in structures:
            self._connection.write(atoms)
