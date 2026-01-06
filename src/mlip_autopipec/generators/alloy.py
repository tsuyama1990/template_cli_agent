# -*- coding: utf-8 -*-
"""
Concrete implementation for creating alloy structures.

This module provides a generator for creating random solid-solution alloy
structures based on a specified composition and lattice type.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase import Atom
from ase.build import bulk

from .base import BaseStructureGenerator

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import SystemConfig


class AlloyGenerator(BaseStructureGenerator):
    """
    A concrete implementation for creating random solid-solution alloy structures.

    This generator creates a supercell of a given lattice, then randomly assigns
    atomic species to the lattice sites according to the specified composition.
    """

    def __init__(self, config: SystemConfig) -> None:
        """
        Initialize the AlloyGenerator.

        Args:
            config: The system configuration Pydantic model.
        """
        self.config = config

    def generate(self) -> list[Atoms]:
        """
        Generate a list of random solid-solution alloy structures.

        The process involves creating a supercell of the primary element and then
        randomly substituting atoms to match the target composition. Each generated
        structure is checked for physical validity.

        Returns:
            A list of ase.Atoms objects.
        """
        structures = []
        for _ in range(self.config.num_structures):
            # Create a primitive cell of the primary element
            primitive = bulk(self.config.elements[0], self.config.lattice, a=4.0, cubic=True)

            # Create a supercell
            num_atoms_primitive = len(primitive)
            target_atoms = 100  # A reasonable number of atoms for a supercell
            scaling_factor = int(np.ceil((target_atoms / num_atoms_primitive) ** (1 / 3)))
            supercell = primitive.repeat((scaling_factor, scaling_factor, scaling_factor))  # type: ignore

            # Get atom counts based on composition
            num_total_atoms = len(supercell)
            atom_counts = {
                element: int(np.round(self.config.composition[element] * num_total_atoms))
                for element in self.config.elements
            }

            # Adjust counts to match total number of atoms due to rounding
            current_total = sum(atom_counts.values())
            if current_total != num_total_atoms:
                diff = num_total_atoms - current_total
                atom_counts[self.config.elements[0]] += diff

            # Create a list of atomic numbers based on counts
            atomic_numbers = []
            for element, count in atom_counts.items():
                atomic_number = Atom(element).number  # type: ignore
                atomic_numbers.extend([atomic_number] * count)

            np.random.shuffle(atomic_numbers)

            # Set the atomic numbers in the supercell
            supercell.set_atomic_numbers(atomic_numbers)

            # Perform a physical validity check
            if self._is_physically_valid(supercell):
                structures.append(supercell)

        return structures

    def _is_physically_valid(self, atoms: Atoms, cutoff: float = 1.0) -> bool:
        """
        Check if the structure is physically valid by checking inter-atomic distances.

        Args:
            atoms: The ase.Atoms object to check.
            cutoff: The minimum allowed distance between atoms in Angstroms.

        Returns:
            True if the structure is valid, False otherwise.
        """
        distances = atoms.get_all_distances(mic=True)  # type: ignore
        min_dist = np.min(distances[np.nonzero(distances)])
        return bool(min_dist > cutoff)
