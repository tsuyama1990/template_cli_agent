from typing import Dict, Any
import ase.db
from ase.atoms import Atoms
from ase.stress import full_3x3_to_voigt_6_stress
from .models import DFTResult

class AseDB:
    def __init__(self, db_path: str):
        self.db_path = db_path

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT calculation results to the database.

        Args:
            atoms: The ase.Atoms object representing the atomic structure.
            result: A DFTResult object containing the calculation outputs.

        Returns:
            The integer ID of the newly created row in the database.
        """
        with ase.db.connect(self.db_path) as db:
            # Prepare key-value pairs from the DFTResult model
            key_value_pairs: Dict[str, Any] = {
                "total_energy_ev": result.total_energy_ev,
                "was_successful": result.was_successful,
            }
            if result.error_message:
                key_value_pairs["error_message"] = result.error_message

            # ASE's SinglePointCalculator is used to store energy, forces, and stress
            # directly with the Atoms object in the database.
            calculator_results = {
                "energy": result.total_energy_ev,
                "forces": result.forces,
                "stress": full_3x3_to_voigt_6_stress(result.stress),
            }

            # Use a SinglePointCalculator to attach results to the Atoms object
            from ase.calculators.singlepoint import SinglePointCalculator
            atoms.calc = SinglePointCalculator(atoms, **calculator_results)

            db_id = db.write(atoms, key_value_pairs=key_value_pairs)
        return db_id

    def get(self, db_id: int) -> tuple[Atoms, Dict[str, Any]]:
        """
        Retrieves an atomic structure and its associated data from the database.

        Args:
            db_id: The ID of the row to retrieve.

        Returns:
            A tuple containing the ase.Atoms object and a dictionary of the key-value pairs.
        """
        with ase.db.connect(self.db_path) as db:
            row = db.get(id=db_id)
            atoms = row.toatoms()
            key_value_pairs = dict(row.key_value_pairs)
        return atoms, key_value_pairs
