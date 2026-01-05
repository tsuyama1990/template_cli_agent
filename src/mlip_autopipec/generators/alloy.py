"""Concrete implementation for creating alloy structures."""
import random
from typing import List

import numpy as np
from ase import Atoms
from ase.build import bulk

from mlip_autopipec.config.models import SystemConfig
from .base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """Generates random solid-solution alloy structures."""

    def __init__(self, config: SystemConfig):
        self.config = config

    def generate(self) -> List[Atoms]:
        """
        Creates a list of random alloy supercells based on the configuration.

        Returns:
            A list of physically valid ase.Atoms objects.
        """
        structures: List[Atoms] = []
        for _ in range(self.config.num_structures):
            # Create a base supercell of the first element.
            # A lattice constant of 4.0 is a reasonable generic value.
            primitive = bulk(self.config.elements[0], self.config.lattice, a=4.0, cubic=True)
            # Using a 3x3x3 supercell to get a reasonable number of atoms.
            supercell = primitive.repeat((3, 3, 3))
            num_atoms = len(supercell)

            # Determine the number of atoms of each element based on composition.
            atom_counts = {el: int(round(num_atoms * frac)) for el, frac in self.config.composition.items()}

            # Adjust counts to match the total number of atoms due to rounding errors.
            current_total = sum(atom_counts.values())
            diff = num_atoms - current_total
            if diff != 0:
                # Add or remove atoms from the element with the largest fractional part
                # to minimize the composition error.
                frac_parts = {el: (num_atoms * frac) % 1 for el, frac in self.config.composition.items()}
                sorted_elements = sorted(frac_parts.keys(), key=lambda el: frac_parts[el], reverse=True)
                for i in range(abs(diff)):
                    el_to_adjust = sorted_elements[i % len(sorted_elements)]
                    atom_counts[el_to_adjust] += np.sign(diff)

            # Create a list of atomic symbols with the correct composition.
            symbols: List[str] = []
            for element, count in atom_counts.items():
                symbols.extend([element] * count)

            # Shuffle the symbols and assign them to the supercell atoms.
            random.shuffle(symbols)
            supercell.set_chemical_symbols(symbols)

            # Ensure the generated structure is physically plausible.
            if self._is_physically_valid(supercell):
                structures.append(supercell)

        return structures

    def _is_physically_valid(self, atoms: Atoms, min_dist: float = 1.0) -> bool:
        """
        Checks if a structure is physically valid by ensuring atoms are not too close.

        Args:
            atoms: The ase.Atoms object to check.
            min_dist: The minimum allowed distance between any two atoms.

        Returns:
            True if the structure is valid, False otherwise.
        """
        distances = atoms.get_all_distances(mic=True)
        # Set the diagonal (self-distances) to infinity to ignore them.
        np.fill_diagonal(distances, np.inf)
        return np.min(distances) > min_dist
