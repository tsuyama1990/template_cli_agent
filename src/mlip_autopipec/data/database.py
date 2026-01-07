from typing import Any

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """
    A wrapper class for the ASE database to handle storing and retrieving atomic data.
    """
    def __init__(self, db_path: str):
        self.db_path = db_path

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT calculation results to the database.

        If the calculation was successful, energy, forces, and stress are attached via
        a SinglePointCalculator. Otherwise, failure metadata is stored in key_value_pairs.

        Args:
            atoms: The ASE Atoms object representing the atomic structure.
            result: A DFTResult object containing the calculation results.
        Returns:
            The ID of the newly written row in the database.
        """
        kvp = {
            "was_successful": result.was_successful,
            "error_message": result.error_message,
        }

        # Clean atoms object of any previous calculator
        atoms.calc = None

        if result.was_successful:
            # ASE's SinglePointCalculator requires stress in the 6-element Voigt form.
            stress_voigt = full_3x3_to_voigt_6_stress(result.stress)

            calculator = SinglePointCalculator(
                atoms,
                energy=result.total_energy_ev,
                forces=result.forces,
                stress=stress_voigt,
            )
            atoms.calc = calculator

        # Filter out None values from metadata
        kvp_filtered = {k: v for k, v in kvp.items() if v not in [None, []]}

        with connect(self.db_path) as db:
            db_id = db.write(atoms, key_value_pairs=kvp_filtered)

        return db_id

    def get(self, db_id: int) -> tuple[Atoms, dict[str, Any]]:
        """
        Retrieves an Atoms object and its metadata from the database.
        Args:
            db_id: The ID of the row to retrieve.
        Returns:
            A tuple containing the ASE Atoms object and a dictionary of the key-value pairs.
            The Atoms object will have a calculator attached if the run was successful.
        """
        try:
            with connect(self.db_path) as db:
                row = db.get(id=db_id)
            if not row:
                # This path is unlikely with db.get but good for robustness
                raise KeyError(f"No entry found with id {db_id}")

            atoms = row.toatoms()
            key_value_pairs = dict(row.key_value_pairs)
            return atoms, key_value_pairs

        except KeyError as e:
            # Re-raise with a more informative message
            raise KeyError(f"No entry found with id {db_id}") from e
