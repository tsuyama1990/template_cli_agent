from typing import Any

from ase import Atoms
from ase.db import connect

from mlip_autopipec.data.models import DFTResult


class AseDB:
    def __init__(self, db_path: str):
        self.db_path = db_path

    def write(self, atoms: Atoms, dft_result: DFTResult) -> int:
        """Writes an Atoms object and its corresponding DFT result to the database."""
        with connect(self.db_path) as db:
            key_value_pairs = dft_result.model_dump()

            # Exclude keys that are stored by the calculator on the Atoms object
            calculator_keys = ["total_energy_ev", "forces", "stress"]
            for key in calculator_keys:
                key_value_pairs.pop(key, None)

            # The ASE db writer can't handle None values for key-value pairs
            key_value_pairs = {k: v for k, v in key_value_pairs.items() if v is not None}

            db_id = db.write(atoms, key_value_pairs=key_value_pairs)
        return db_id

    def get(self, db_id: int) -> tuple[Atoms, dict[str, Any]] | None:
        """Retrieves an Atoms object and its key-value pairs from the database."""
        with connect(self.db_path) as db:
            try:
                row = db.get(id=db_id)
                atoms = row.toatoms()
                key_value_pairs = row.key_value_pairs
                return atoms, key_value_pairs
            except KeyError:
                return None
