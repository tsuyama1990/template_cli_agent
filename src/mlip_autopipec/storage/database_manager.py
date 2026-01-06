# src/mlip_autopipec/storage/database_manager.py
"""
Manages all interactions with the ASE (Atomic Simulation Environment) database,
which is used to store the final, sampled atomic structures.
"""

from __future__ import annotations

from collections.abc import Generator
from contextlib import contextmanager
from typing import TYPE_CHECKING

from ase.db import connect
from ase.db.core import Database

if TYPE_CHECKING:
    from ase import Atoms


class DatabaseManager:
    """Handles writing atomic structures to an ASE database."""

    @contextmanager
    def connect(self, path: str) -> Generator[Database, None, None]:
        """
        Provides a context manager for connecting to the ASE database.
        Ensures the connection is properly closed.

        Args:
            path: The file path to the database.

        Yields:
            The ASE database connection object.
        """
        # The ase.db.connect object doesn't need to be explicitly closed.
        # It handles its connection internally.
        connection = connect(path)  # type: ignore[no-untyped-call]
        yield connection

    def write_structures(self, connection: Database, structures: list[Atoms]) -> None:
        """
        Writes a list of ASE Atoms objects to the given database connection.

        Args:
            connection: The active ASE database connection.
            structures: A list of structures to write.
        """
        for atoms in structures:
            connection.write(atoms)
