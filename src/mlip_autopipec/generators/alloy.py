"""Concrete implementation for creating alloy structures."""

import random

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """A generator for creating alloy structures."""

    def __init__(self, config: SystemConfig) -> None:
        self.config = config

    def generate(self) -> list[Atoms]:
        """
        Generates a list of random solid-solution alloy structures.

        Returns:
            A list of ASE Atoms objects.
        """
        structures = []
        for _ in range(self.config.num_structures):
            primitive_cell = bulk(self.config.elements[0], self.config.lattice, a=3.6)
            supercell_size = self._get_supercell_size(len(primitive_cell))
            supercell = primitive_cell.repeat(supercell_size)

            num_atoms = len(supercell)
            symbols = []
            for element, fraction in self.config.composition.items():
                num_element_atoms = round(fraction * num_atoms)
                symbols.extend([element] * num_element_atoms)

            # Adjust for rounding errors
            while len(symbols) < num_atoms:
                symbols.append(random.choice(self.config.elements))  # noqa: S311

            random.shuffle(symbols)
            supercell.set_chemical_symbols(symbols)

            if self._is_physically_valid(supercell):
                structures.append(supercell)

        return structures

    def _get_supercell_size(self, num_primitive_atoms: int) -> int:
        """Determines the supercell size to get at least 64 atoms."""
        return int(np.ceil((64 / num_primitive_atoms) ** (1 / 3)))

    def _is_physically_valid(self, atoms: Atoms) -> bool:
        """Checks if the structure is physically valid."""
        distances = atoms.get_all_distances(mic=True)
        min_dist = np.min(distances[np.nonzero(distances)])
        return min_dist > 1.0
