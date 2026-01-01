"""
This module implements the structure generator for metallic alloys.

It uses a randomized supercell approach to generate a diverse set of initial
structures that approximate a disordered alloy, serving as a robust starting
point for the MLIP workflow.
"""

import logging

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.build.supercells import make_supercell
from pymatgen.core import Composition

from ..constants import ALLOY_STRAIN_MAGNITUDE, ALLOY_TARGET_ATOMS
from .base_generator import BaseStructureGenerator

logger = logging.getLogger(__name__)


class AlloyGenerator(BaseStructureGenerator):
    """
    Generates initial structures for metallic alloys.

    This generator creates a randomized, disordered supercell to represent the
    alloy and then applies various strains to generate a diverse dataset.
    """

    def generate(self) -> list[Atoms]:
        """
        Main method to generate a list of alloy structures.

        Returns:
            A list of ase.Atoms objects, including one base structure and
            several strained variations.
        """
        base_atoms = self._create_randomized_supercell()
        strained_atoms_list = self._apply_strains(base_atoms)
        return [base_atoms] + strained_atoms_list

    def _create_randomized_supercell(self) -> Atoms:
        """
        Creates a single, large, randomized supercell for the alloy.

        This method builds a simple cubic primitive cell, scales it to a
        supercell of a target size, and then randomly assigns atomic species
        to the lattice sites according to the system's composition.

        Returns:
            An ase.Atoms object representing the disordered supercell.
        """
        primitive_cell = bulk(self.system_info.elements[0], "sc", a=3.5, cubic=True)

        scaling_factor = int(np.ceil((ALLOY_TARGET_ATOMS / len(primitive_cell)) ** (1 / 3)))
        scaling_matrix = np.diag([scaling_factor] * 3)
        supercell = make_supercell(primitive_cell, scaling_matrix)

        composition = Composition(self.system_info.composition).fractional_composition
        total_atoms = len(supercell)
        atomic_numbers = self._get_atomic_numbers_from_composition(
            composition, total_atoms
        )

        np.random.shuffle(atomic_numbers)
        supercell.set_atomic_numbers(atomic_numbers)

        logger.info(f"Created a randomized supercell with formula: {supercell.get_chemical_formula()}")
        return supercell

    def _get_atomic_numbers_from_composition(
        self, composition: Composition, total_atoms: int
    ) -> list[int]:
        """
        Determines the exact number of atoms of each element.
        """
        atomic_numbers = []
        for element, fraction in composition.items():
            count = round(fraction * total_atoms)
            atomic_numbers.extend([element.Z] * int(count))

        # Adjust counts to match the total number of atoms due to rounding
        while len(atomic_numbers) < total_atoms:
            atomic_numbers.append(composition.elements[-1].Z)

        return atomic_numbers[:total_atoms]

    def _apply_strains(self, base_atoms: Atoms) -> list[Atoms]:
        """
        Applies random volumetric and shear strains to a base structure.
        """
        strained_structures = []
        num_strains = self.config.simulation.initial_structures_to_generate - 1
        for _ in range(num_strains):
            strained_atoms = base_atoms.copy()
            strain = np.random.uniform(
                -ALLOY_STRAIN_MAGNITUDE, ALLOY_STRAIN_MAGNITUDE, (3, 3)
            )
            strain = (strain + strain.T) / 2  # Symmetrize
            new_cell = strained_atoms.get_cell() @ (np.identity(3) + strain)
            strained_atoms.set_cell(new_cell, scale_atoms=True)
            strained_structures.append(strained_atoms)

        logger.info(f"Applied {num_strains} random strains to the base structure.")
        return strained_structures
