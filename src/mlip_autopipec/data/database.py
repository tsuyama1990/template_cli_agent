"""Interface to the ASE database for storing atomic structures and results."""
from pathlib import Path
from typing import Union, Tuple

import ase.db
import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """A wrapper around the ASE database for project-specific operations."""

    def __init__(self, db_path: Union[str, Path]):
        """
        Initializes the database connection.

        Args:
            db_path: Path to the ASE database file.
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._connection = ase.db.connect(self.db_path)

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFT result to the database.

        If the calculation was successful, the results (energy, forces, stress)
        are attached to the Atoms object via a SinglePointCalculator.
        Metadata is stored in the key-value pairs.

        Args:
            atoms: The ase.Atoms object.
            result: The DFTResult Pydantic model.

        Returns:
            The unique ID of the new row in the database.
        """
        metadata = {
            "was_successful": result.was_successful,
            "error_message": result.error_message,
        }
        metadata = {k: v for k, v in metadata.items() if v is not None}

        if result.was_successful:
            # ASE's SinglePointCalculator expects stress in the 6-element Voigt form.
            stress_voigt = None
            if result.stress is not None:
                stress_voigt = full_3x3_to_voigt_6_stress(np.array(result.stress))

            atoms.calc = SinglePointCalculator(
                atoms,
                energy=result.total_energy_ev,
                forces=result.forces,
                stress=stress_voigt,
            )
        elif hasattr(atoms, 'calc'):
            atoms.calc = None

        db_id = self._connection.write(atoms, key_value_pairs=metadata)
        return db_id

    def get(self, db_id: int) -> Tuple[Atoms, DFTResult]:
        """
        Retrieves a record from the database by its ID.

        Args:
            db_id: The unique ID of the row to retrieve.

        Returns:
            A tuple containing the ase.Atoms object and the reconstructed
            DFTResult Pydantic model.
        """
        row = self._connection.get(db_id)
        atoms = row.toatoms()
        metadata = row.key_value_pairs

        if metadata.get("was_successful"):
            # Use the high-level ASE APIs to retrieve results.
            # This correctly handles conversions (e.g., Voigt stress to 3x3).
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces() if atoms.calc.results.get('forces') is not None else None
            stress = atoms.get_stress(voigt=False) if atoms.calc.results.get('stress') is not None else None

            dft_result = DFTResult(
                total_energy_ev=energy,
                forces=forces.tolist() if forces is not None else None,
                stress=stress.tolist() if stress is not None else None,
                was_successful=True,
            )
        else:
            dft_result = DFTResult(
                was_successful=False,
                error_message=metadata.get("error_message"),
            )

        return atoms, dft_result
