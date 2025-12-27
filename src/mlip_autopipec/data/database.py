from ase.db import connect
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
from typing import Dict, Any, Tuple

from mlip_autopipec.data.models import DFTResult

class AseDB:
    def __init__(self, db_path: str):
        self.db_path = db_path

    def write(self, atoms: Atoms, dft_result: DFTResult) -> int:
        """Writes an ASE Atoms object and its corresponding DFT result to the database."""
        # ASE DB reserves keys like 'energy', 'forces', 'stress'.
        # These must be attached via a calculator.
        if dft_result.was_successful:
            calc = SinglePointCalculator(
                atoms,
                energy=dft_result.total_energy_ev,
                forces=np.array(dft_result.forces),
                stress=np.array(dft_result.stress)
            )
            atoms.calc = calc

        # Other metadata can be stored in key_value_pairs
        key_value_pairs = {
            "was_successful": dft_result.was_successful,
            "error_message": dft_result.error_message,
        }
        # Filter out None values
        key_value_pairs = {k: v for k, v in key_value_pairs.items() if v is not None}

        with connect(self.db_path) as db:
            db_id = db.write(atoms, key_value_pairs=key_value_pairs)

        return db_id

    def get(self, db_id: int) -> Tuple[Atoms, Dict[str, Any]]:
        """Retrieves an Atoms object and its key-value pairs from the database."""
        with connect(self.db_path) as db:
            row = db.get(id=db_id)
            atoms = row.toatoms()
            key_value_pairs = row.key_value_pairs
        return atoms, key_value_pairs
