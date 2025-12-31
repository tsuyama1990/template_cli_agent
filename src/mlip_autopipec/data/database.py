"""Database wrapper for ASE DB."""

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """A wrapper class for the ASE database to abstract database operations."""

    def __init__(self, db_path: str):
        """
        Initializes the AseDB wrapper and connects to the database.

        Args:
            db_path: Path to the ASE database file.
        """
        self.db_path = db_path
        self.connection = connect(self.db_path)

    def add_atoms(self, atoms_list: list[Atoms]):
        """
        Adds a list of ASE Atoms objects to the database.

        Args:
            atoms_list: A list of Atoms objects to add.
        """
        for atoms in atoms_list:
            self.connection.write(atoms, state="unlabelled")

    def get_atoms_by_state(self, state: str) -> list[Atoms]:
        """
        Retrieves a list of Atoms objects from the database by their state.

        Args:
            state: The state to filter by (e.g., 'unlabelled', 'labelled').

        Returns:
            A list of Atoms objects matching the state.
        """
        atoms_list = []
        for row in self.connection.select(f"state={state}"):
            atoms = row.toatoms()
            atoms.info["db_id"] = row.id  # Attach db id for later reference
            atoms_list.append(atoms)
        return atoms_list

    def update_with_dft_results(self, atoms_id: int, results: DFTResult):
        """
        Updates a database entry with the results of a DFT calculation.

        This method retrieves the relevant Atoms object, attaches a new
        SinglePointCalculator with the provided results, and then updates
        the database entry.

        Args:
            atoms_id: The unique ID of the atoms entry in the database.
            results: A DFTResult object containing the calculated data.
        """
        row = self.connection.get(id=atoms_id)
        atoms = row.toatoms()

        stress_voigt = results.stress
        if results.stress.shape == (3, 3):
            stress_voigt = full_3x3_to_voigt_6_stress(results.stress)

        calc = SinglePointCalculator(
            atoms=atoms,
            energy=results.energy,
            forces=results.forces,
            stress=stress_voigt,
        )
        atoms.calc = calc

        self.connection.update(id=atoms_id, atoms=atoms, state="labelled")
