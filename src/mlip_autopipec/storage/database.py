# -*- coding: utf-8 -*-
"""
AseDBWrapper for abstracting database interactions.

This module provides a context manager for safe and convenient interaction
with ASE databases.
"""
from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import TYPE_CHECKING, Optional, Type

from ase.db import connect

if TYPE_CHECKING:
    from ase import Atoms
    from ase.db.core import Database


class AseDBWrapper:
    """
    A context manager to handle the creation of and connection to the output SQLite database.

    This wrapper simplifies database interactions by managing the connection
    lifecycle, ensuring that the database is correctly connected before use and
    that resources are properly handled.

    Usage:
        with AseDBWrapper(Path("my_database.db")) as db:
            db.write_structures(list_of_atoms)
    """

    def __init__(self, db_path: Path) -> None:
        """
        Initialize the AseDBWrapper.

        Args:
            db_path: The path to the database file.
        """
        self.db_path = db_path
        self._db: Optional[Database] = None

    def __enter__(self) -> "AseDBWrapper":
        """
        Enter the context manager and connect to the database.

        Returns:
            The instance of the AseDBWrapper.
        """
        self._db = connect(self.db_path)  # type: ignore
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        """
        Exit the context manager and ensure the connection is no longer accessible.

        This method handles the cleanup of the database connection.
        """
        if self._db:
            # The db object from ase.db.connect doesn't have a close method,
            # it's managed by Python's garbage collection.
            # We just dereference it to signify disconnection.
            self._db = None

    def write_structures(self, structures: list[Atoms]) -> None:
        """
        Write a list of Atoms objects to the database.

        Args:
            structures: A list of ase.Atoms objects to write.

        Raises:
            RuntimeError: If the database is not connected (i.e., not used within a 'with' block).
        """
        if self._db is None:
            msg = "Database not connected. Use this class as a context manager."
            raise RuntimeError(msg)
        for atoms in structures:
            self._db.write(atoms)
