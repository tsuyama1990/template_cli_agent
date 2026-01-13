from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.db import connect

if TYPE_CHECKING:
    from ase import Atoms


class AseDBWrapper:
    """A wrapper class for the ASE database to handle database connections."""

    def __init__(self, db_path: Path):
        self.db_path = db_path

    def write_structures(self, structures: list[Atoms]) -> None:
        """Writes a list of ASE Atoms objects to the database.

        Args:
            structures: A list of ASE Atoms objects to write to the database.
        """
        with connect(self.db_path) as db:
            for atoms in structures:
                db.write(atoms)
