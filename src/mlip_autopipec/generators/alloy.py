from __future__ import annotations

import random
from typing import TYPE_CHECKING

import numpy as np
from ase.build import bulk
from ase.data import atomic_numbers

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.base import BaseStructureGenerator

if TYPE_CHECKING:
    from ase import Atoms


class AlloyGenerator(BaseStructureGenerator):
    """A generator for creating alloy structures."""

    def __init__(self, config: SystemConfig):
        self.config = config

    def generate(self) -> list[Atoms]:
        """Generates a list of random alloy structures."""
        structures = []
        for _ in range(self.config.num_structures):
            # Create a primitive cell
            primitive_cell = bulk(
                self.config.elements[0], self.config.lattice, a=3.6, cubic=True
            )

            # Create a supercell
            num_atoms = len(primitive_cell)
            supercell_size = int(np.ceil((64 / num_atoms) ** (1 / 3)))
            supercell = primitive_cell * (supercell_size, supercell_size, supercell_size)

            # Get the total number of atoms in the supercell
            total_atoms = len(supercell)

            # Determine the number of atoms of each element
            num_atoms_per_element = {}
            remaining_atoms = total_atoms
            for element, fraction in self.config.composition.items():
                num_atoms = round(fraction * total_atoms)
                num_atoms_per_element[element] = num_atoms
                remaining_atoms -= num_atoms

            # Distribute remaining atoms due to rounding
            if remaining_atoms > 0:
                elements = list(self.config.composition.keys())
                for _ in range(remaining_atoms):
                    element_to_add = random.choice(elements)
                    num_atoms_per_element[element_to_add] += 1

            # Create a list of atomic numbers for the supercell
            atomic_numbers_list = []
            for element, num_atoms in num_atoms_per_element.items():
                atomic_numbers_list.extend(
                    [self._get_atomic_number(element)] * num_atoms
                )

            # Shuffle the atomic numbers and set them in the supercell
            random.shuffle(atomic_numbers_list)
            supercell.set_atomic_numbers(atomic_numbers_list)

            # Physical validation to prevent overlapping atoms
            distances = supercell.get_all_distances(mic=True)
            min_distance = np.min(distances[np.nonzero(distances)])
            if min_distance < 1.0:
                # Retry if atoms are too close, this is a simplified approach
                # A more robust solution might involve small random displacements
                continue

            structures.append(supercell)
        return structures

    def _get_atomic_number(self, element: str) -> int:
        """Helper to get atomic number."""
        return atomic_numbers[element]
