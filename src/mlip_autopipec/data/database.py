from pathlib import Path
from typing import Any

import ase.db
from ase.atoms import Atoms

from .models import DFTResult


class AseDB:
    """A wrapper class for the ASE database."""

    def __init__(self, db_path: Path):
        """
        Initializes the AseDB wrapper.

        Args:
            db_path: The path to the ASE database file.
        """
        self._db_path = db_path
        # Explicitly set type='db' to avoid issues with ASE's file type detection
        self._connection = ase.db.connect(self._db_path, type="db")

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT result to the database.

        Args:
            atoms: The ASE Atoms object.
            result: The DFTResult Pydantic model.

        Returns:
            The ID of the new row in the database.
        """
        key_value_pairs = result.dict()
        # Filter out None values to prevent ASE ValueError
        key_value_pairs = {k: v for k, v in key_value_pairs.items() if v is not None}
        return self._connection.write(atoms, key_value_pairs=key_value_pairs)

    def get(self, selection: Any) -> tuple[Atoms | None, dict[str, Any] | None]:
        """
        Retrieves a row from the database.

        Args:
            selection: The selection criteria for the row.

        Returns:
            A tuple containing the Atoms object and the key-value pairs, or (None, None) if not found.
        """
        try:
            row = self._connection.get(selection)
            if row:
                return row.toatoms(), row.key_value_pairs
            return None, None
        except (KeyError, IndexError):
            return None, None

    @property
    def connection(self):
        return self._connection
