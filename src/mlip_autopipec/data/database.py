from ase.db import connect
from ase.atoms import Atoms
from typing import Dict, Any, Optional

class AseDB:
    def __init__(self, db_path: str):
        self.db_path = db_path
        self._validate_db_path()

    def _validate_db_path(self):
        if not self.db_path.endswith('.db'):
            raise ValueError("Database path must end with '.db'")

    def write(self, atoms: Atoms, metadata: Optional[Dict[str, Any]] = None) -> int:
        """Writes the atoms object and its metadata to the database."""
        with connect(self.db_path) as db:
            # Filter out None values to prevent ASE db errors
            filtered_metadata = {k: v for k, v in metadata.items() if v is not None} if metadata else {}
            db_id = db.write(atoms, key_value_pairs=filtered_metadata)
        return db_id

    def read(self, db_id: int) -> Dict[str, Any]:
        """Reads a row from the database by its ID."""
        with connect(self.db_path) as db:
            row = db.get(db_id)
        return row.key_value_pairs
