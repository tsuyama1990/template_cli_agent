# src/mlip_autopipec/storage/database_manager.py
"""Manages all interactions with the ASE database."""

from typing import List
from ase import Atoms
from ase.db import connect


class DatabaseManager:
    """Provides a clean interface for writing to an ASE database."""

    def __init__(self, db_path: str) -> None:
        self.db_path = db_path

    def write_structures(self, structures: list[Atoms]) -> None:
        """
        Writes a list of ASE Atoms objects to the database.

        Energy, forces, and stress are read from the calculator attached
        to each Atoms object.

        Args:
            structures: The list of Atoms objects to write.
        """
        with connect(self.db_path) as db:  # type: ignore[no-untyped-call]
            for atoms in structures:
                db.write(atoms)  # type: ignore[no-untyped-call]
