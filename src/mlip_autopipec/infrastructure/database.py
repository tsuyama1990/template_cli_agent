"""Concrete implementation of the DatabasePort using ASE DB."""

from typing import List
from ase.db import connect
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.domain.models import DFTResult
from mlip_autopipec.domain.ports import DatabasePort


class AseDB(DatabasePort):
    """A wrapper class for the ASE database that implements the DatabasePort."""

    def __init__(self, db_path: str):
        """
        Initializes the AseDB wrapper.

        Args:
            db_path: The path to the database file.
        """
        self.db_path = db_path
        self.connection = connect(self.db_path)

    def add_atoms(self, atoms_list: List[Atoms]):
        for atoms in atoms_list:
            self.connection.write(atoms, state='unlabelled')

    def get_atoms_by_state(self, state: str) -> List[Atoms]:
        atoms_list = []
        for row in self.connection.select(f'state={state}'):
            atoms = row.toatoms()
            atoms.info['db_id'] = row.id
            atoms_list.append(atoms)
        return atoms_list

    def update_with_dft_results(self, atoms_id: int, results: DFTResult):
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

        self.connection.update(id=atoms_id, atoms=atoms, state='labelled')
