"""Interface to the ASE database for storing and retrieving calculation results."""
from typing import Dict, Any, Optional
from typing import Dict, Any, Optional
from ase.db import connect
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress
from .models import DFTResult
import numpy as np

class AseDB:
    """A wrapper class for the ASE database."""

    def __init__(self, db_path: str):
        """
        Initializes the AseDB wrapper.

        Args:
            db_path: Path to the ASE database file.
        """
        self.db_path = db_path

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT result to the database.

        Args:
            atoms: The ASE Atoms object.
            result: The DFTResult object.

        Returns:
            The ID of the new database entry.
        """
        with connect(self.db_path) as db:
            # Attach results to atoms object via a calculator
            # to be compliant with ASE db schema.
            forces_np = np.array(result.forces)
            if forces_np.shape != (len(atoms), 3):
                forces_np = np.zeros((len(atoms), 3))

            stress_np = np.array(result.stress)
            if stress_np.shape != (3, 3):
                stress_np = np.zeros((3, 3))
            stress_voigt = full_3x3_to_voigt_6_stress(stress_np)

            atoms.calc = SinglePointCalculator(
                atoms=atoms,
                energy=result.total_energy_ev,
                forces=forces_np,
                stress=stress_voigt,
            )

            key_value_pairs = {
                "was_successful": result.was_successful,
                "error_message": result.error_message,
            }
            # Filter out None values to prevent db errors
            key_value_pairs = {k: v for k, v in key_value_pairs.items() if v is not None}

            db_id = db.write(atoms, key_value_pairs=key_value_pairs)
        return db_id

    def get(self, db_id: int) -> Optional[Dict[str, Any]]:
        """
        Retrieves a record from the database by its ID.

        Args:
            db_id: The ID of the database entry.

        Returns:
            A dictionary containing the retrieved data, or None if not found.
        """
        with connect(self.db_path) as db:
            try:
                row = db.get(id=db_id)
            except KeyError:
                return None

        if not row:
            return None

        atoms = row.toatoms()
        return {
            "atoms": atoms,
            "total_energy_ev": atoms.get_potential_energy(),
            "forces": atoms.get_forces(),
            "stress": atoms.get_stress(voigt=False),
            "was_successful": row.key_value_pairs.get("was_successful", False),
            "error_message": row.key_value_pairs.get("error_message"),
        }
