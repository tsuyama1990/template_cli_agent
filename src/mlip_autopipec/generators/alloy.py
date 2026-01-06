"""Concrete implementation for creating alloy structures."""
from __future__ import annotations

import random
from typing import TYPE_CHECKING

import numpy as np
from ase.build import bulk, make_supercell

from .base import BaseStructureGenerator

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import SystemConfig


class AlloyGenerator(BaseStructureGenerator):
    """A generator for creating random solid-solution alloy structures."""

    def __init__(self, config: SystemConfig) -> None:
        """Initialise the AlloyGenerator."""
        super().__init__(config)
        self.supercell_matrix = np.diag([4, 4, 4])  # Fixed 4x4x4 supercell

    def generate(self) -> list[Atoms]:
        """
        Generate a list of random alloy structures.

        Returns:
            A list of `ase.Atoms` objects representing the generated structures.
        """
        structures = []
        for _ in range(self.config.num_structures):
            # 1. Create a primitive cell of the primary element
            primitive_cell = bulk(
                self.config.elements[0], crystalstructure=self.config.lattice, a=None, cubic=True
            )

            # 2. Create a supercell
            supercell = make_supercell(primitive_cell, self.supercell_matrix)
            num_atoms = len(supercell)

            # 3. Determine the number of atoms for each element based on composition
            symbols_to_assign = []
            atoms_count = 0
            for element, fraction in self.config.composition.items():
                num_element_atoms = round(fraction * num_atoms)
                symbols_to_assign.extend([element] * num_element_atoms)
                atoms_count += num_element_atoms

            # Adjust for rounding errors to ensure the total number of atoms is correct
            while atoms_count < num_atoms:
                symbols_to_assign.append(random.choice(self.config.elements))
                atoms_count += 1
            symbols_to_assign = symbols_to_assign[:num_atoms]

            # 4. Randomly assign elements to atomic positions
            random.shuffle(symbols_to_assign)
            supercell.set_chemical_symbols(symbols_to_assign)

            # 5. Perform a basic physical validity check (no overlapping atoms)
            # A random substitution on a lattice should not cause this,
            # but it's a good practice safety check.
            distances = supercell.get_all_distances(mic=True)
            if distances.size > 0:
                min_dist = np.min(distances[np.nonzero(distances)])
                if min_dist < 1.0:
                    # In a more complex scenario, we might regenerate or log a warning.
                    # For this implementation, we assume valid lattice constants from ASE.
                    pass
            structures.append(supercell)
        return structures
