# Description: Interface to the ASE database for storing and retrieving atomic structures and calculation results.
from typing import Any

import ase.db
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """A wrapper class for the ASE database to handle I/O operations."""

    def __init__(self, path: str):
        """
        Initializes the database connection.

        Args:
            path: The file path to the ASE database.
        """
        self.db = ase.db.connect(path)

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an atomic structure and its corresponding DFT results to the database.

        If the calculation was successful, the energy, forces, and stress are attached
        to the Atoms object via a SinglePointCalculator. Other metadata is stored
        in the key-value pairs.

        Args:
            atoms: The ase.Atoms object representing the structure.
            result: A DFTResult object containing the calculation results.

        Returns:
            The integer ID of the newly created row in the database.
        """
        kvp = {
            "was_successful": result.was_successful,
            "error_message": result.error_message,
        }

        if result.was_successful:
            voigt_stress = (
                full_3x3_to_voigt_6_stress(result.stress) if result.stress else None
            )
            calc = SinglePointCalculator(
                atoms,
                energy=result.total_energy_ev,
                forces=result.forces,
                stress=voigt_stress,
            )
            atoms.calc = calc

        # Filter out None values as ase.db raises a ValueError for them.
        kvp_filtered = {k: v for k, v in kvp.items() if v is not None}
        db_id = self.db.write(atoms, key_value_pairs=kvp_filtered)
        return db_id

    def get(self, db_id: int) -> tuple[Atoms, dict[str, Any]]:
        """
        Retrieves an atomic structure and its metadata from the database.

        Args:
            db_id: The integer ID of the row to retrieve.

        Returns:
            A tuple containing the ase.Atoms object and a dictionary of the
            key-value pairs.
        """
        row = self.db.get(db_id)
        atoms = row.toatoms()
        # The calculator results are automatically attached to the Atoms object.
        return atoms, row.key_value_pairs
