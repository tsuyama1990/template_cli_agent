"""AseDBWrapper for abstracting database interactions."""
from pathlib import Path
from typing import List

from ase import Atoms
from ase.db import connect


class AseDBWrapper:
    """A wrapper for the ASE database."""

    def __init__(self, db_path: Path):
        """Initialize the AseDBWrapper."""
        self.db_path = db_path

    def write_structures(self, structures: List[Atoms]):
        """Write structures to the database."""
        with connect(self.db_path) as db:
            for atoms in structures:
                db.write(atoms)
