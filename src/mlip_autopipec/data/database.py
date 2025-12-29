# ruff: noqa: D101, D102, D103, D104, D105, D107
from pathlib import Path
from typing import Any

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """A wrapper class for interacting with an ASE (Atomic Simulation Environment)
    SQLite database to store and retrieve atomic structures and calculation results.
    """

    def __init__(self, db_path: Path | str):
        """Initializes the AseDB wrapper.

        Args:
            db_path: The path to the SQLite database file.
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT result to the database.

        Args:
            atoms: The ASE Atoms object to write.
            result: The DFTResult Pydantic model containing calculation results.

        Returns:
            The ID of the newly written database row.
        """
        atoms_to_write = atoms.copy()

        if result.was_successful:
            # Convert 3x3 stress matrix to 6-element Voigt vector for ASE
            voigt_stress = full_3x3_to_voigt_6_stress(result.stress)

            # Attach a SinglePointCalculator to store results on the Atoms object
            calculator = SinglePointCalculator(
                atoms_to_write,
                energy=result.total_energy_ev,
                forces=result.forces,
                stress=voigt_stress,
            )
            atoms_to_write.calc = calculator

        # Key-value pairs should only contain metadata, not calculator results
        key_value_pairs: dict[str, Any] = {
            "was_successful": result.was_successful,
        }
        if result.error_message:
            key_value_pairs["error_message"] = result.error_message

        with connect(self.db_path) as db:
            db_id = db.write(atoms_to_write, key_value_pairs=key_value_pairs)

        return db_id

    def get(self, db_id: int) -> tuple[Atoms, dict[str, Any]]:
        """
        Retrieves an Atoms object and its key-value pairs from the database.

        Args:
            db_id: The ID of the database row to retrieve.

        Returns:
            A tuple containing the ASE Atoms object and its key-value pairs.
        """
        with connect(self.db_path) as db:
            row = db.get(db_id)
            atoms = row.toatoms()
            key_value_pairs = row.key_value_pairs
        return atoms, key_value_pairs
