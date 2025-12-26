import sqlite3
from typing import Optional

from ase.db import connect
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTResult

class AseDB:
    """A wrapper around the ASE database for storing and retrieving calculation results."""

    def __init__(self, db_path: str):
        """
        Initializes the database connection.

        Args:
            db_path: Path to the SQLite database file.
        """
        self.db_path = db_path
        try:
            self.connection = connect(db_path)
        except sqlite3.OperationalError as e:
            raise ConnectionError(f"Could not connect to database at {db_path}. Ensure the directory exists.") from e


    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes a new entry to the database.

        Args:
            atoms: The atomic structure.
            result: The DFT calculation result object.

        Returns:
            The ID of the newly created database row.
        """
        db_id = self.connection.write(atoms, data=result.model_dump())
        return db_id

    def get_row(self, id: int) -> Optional[dict]:
        """
        Retrieves a single row from the database by its ID.

        Args:
            id: The unique ID of the database row.

        Returns:
            A dictionary containing the row data, or None if not found.
        """
        try:
            row = self.connection.get(id=id)
            return row.data
        except KeyError:
            return None

    def get_atoms(self, id: int) -> Optional[Atoms]:
        """
        Retrieves the Atoms object for a given ID.

        Args:
            id: The unique ID of the database row.

        Returns:
            The ase.Atoms object, or None if not found.
        """
        try:
            row = self.connection.get(id=id)
            return row.toatoms()
        except KeyError:
            return None
