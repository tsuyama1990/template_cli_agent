from pathlib import Path

from ase import Atoms
from ase.db import connect


class DatabaseManager:
    """
    Manages all interactions with the ASE database for storing final structures.
    """

    def __init__(self, db_path: Path):
        self.db_path = db_path
        self._connection = None

    def __enter__(self):
        """Connect to the database when entering a context."""
        self._connection = connect(self.db_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Ensure the connection is closed when exiting the context."""
        self._connection = None # ASE's connect() handles closing

    def write_structures(self, structures: list[Atoms]):
        """
        Writes a list of ASE Atoms objects to the database.

        Each structure's energy, forces, and other calculator results are
        automatically included by ASE.

        Args:
            structures: A list of structures to be saved.
        """
        if not self._connection:
            raise ConnectionError("Database is not connected. Use 'with' statement.")

        for atoms in structures:
            self._connection.write(atoms)
