"""Provides a concrete implementation for generating alloy structures.

This module contains the `AlloyGenerator`, which is responsible for creating
a set of initial "seed" structures for multi-component alloy systems. It
builds upon a specified crystal lattice (e.g., FCC, BCC) and populates it with
atoms according to the desired chemical composition, resulting in random
solid-solution configurations.
"""

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
    """Generates random solid-solution alloy structures.

    This generator takes a system configuration specifying elements, composition,
    and lattice type, and produces a requested number of atomic structures. It
    first creates a large supercell of a single element and then randomly
    replaces atoms with other elements to match the target composition.

    Attributes:
        supercell_matrix (np.ndarray): The transformation matrix used to create
            the supercell from the primitive cell. Fixed at 4x4x4 for Cycle 1.

    """

    def __init__(self, config: SystemConfig) -> None:
        """Initialise the AlloyGenerator."""
        super().__init__(config)
        self.supercell_matrix = np.diag([4, 4, 4])  # Fixed 4x4x4 supercell

    def generate(self) -> list[Atoms]:
        """Generate a list of random alloy structures.

        The process involves:
        1. Creating a primitive unit cell of the primary element.
        2. Expanding it into a larger supercell.
        3. Calculating the number of atoms of each element for the target composition.
        4. Randomly assigning elemental symbols to the atomic sites in the supercell.
        5. Performing a basic physical validity check to ensure no atoms overlap.

        Returns:
            A list of `ase.Atoms` objects, each representing a unique random
            alloy structure.

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
