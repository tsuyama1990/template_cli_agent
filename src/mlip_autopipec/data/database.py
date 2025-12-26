from ase.atoms import Atoms
from ase.db import connect
from pathlib import Path
from typing import Union

from mlip_autopipec.data.models import DFTResult

class AseDB:
    """A wrapper class for ASE's database functionality."""

    def __init__(self, db_path: Union[str, Path]):
        """
        Initializes the database connection.

        Args:
            db_path: The path to the SQLite database file.
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._db = connect(self.db_path)

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its associated DFT result to the database.

        Args:
            atoms: The ASE Atoms object.
            result: The DFTResult Pydantic model.

        Returns:
            The unique ID of the new row in the database.
        """
        key_value_pairs = {"was_successful": result.was_successful}
        if result.was_successful:
            key_value_pairs["total_energy_ev"] = result.total_energy_ev

        db_id = self._db.write(
            atoms,
            key_value_pairs=key_value_pairs,
            data=result.model_dump()
        )
        return db_id

    def get(self, db_id: int):
        """
        Retrieves a row from the database by its ID.

        Args:
            db_id: The unique ID of the row.

        Returns:
            An ASE database row object.
        """
        return self._db.get(id=db_id)
