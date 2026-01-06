from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from ase.db import connect

if TYPE_CHECKING:
    from ase import Atoms

logger = logging.getLogger(__name__)


class AseDBWrapper:
    """A wrapper for the ASE database to handle storage operations."""

    def __init__(self, db_path: Path) -> None:
        """
        Initialise the AseDBWrapper.

        Args:
            db_path: The path to the ASE database file.
        """
        self.db_path = db_path
        self.connection = None

    def __enter__(self) -> AseDBWrapper:
        """Enter the context manager and connect to the database."""
        self.connection = connect(self.db_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit the context manager."""
        # The ASE db connection object handles its own closing implicitly.
        pass

    def write_structures(self, structures: list[Atoms]) -> None:
        """
        Write a list of Atoms objects to the database.

        Args:
            structures: A list of ASE Atoms objects to write.
        """
        if self.connection is None:
            raise ConnectionError("Database connection is not open. Use 'with AseDBWrapper(...) as db:'")
        for atoms in structures:
            self.connection.write(atoms)
        logger.info(f"Wrote {len(structures)} structures to {self.db_path}")
