"""Concrete implementation for creating alloy structures."""
import random
from typing import List

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """Generates random solid-solution alloy structures."""

    def __init__(self, config: SystemConfig):
        """
        Initialize the generator with system configuration.

        Args:
            config: The SystemConfig Pydantic model.
        """
        self.config = config

    def generate(self) -> List[Atoms]:
        """
        Generate a list of random alloy structures.

        Returns:
            A list of ase.Atoms objects.
        """
        structures = []
        for _ in range(self.config.num_structures):
            # Use the first element as a template for the bulk structure
            template_element = self.config.elements[0]
            atoms = bulk(template_element, self.config.lattice, a=4.0, cubic=True)

            # Create a supercell to have enough atoms for the composition
            num_atoms_primitive = len(atoms)
            # A simple heuristic to get a decent-sized supercell
            supercell_size = int(np.ceil((50 / num_atoms_primitive) ** (1 / 3)))
            atoms = atoms.repeat((supercell_size, supercell_size, supercell_size))

            # Get the total number of atoms and assign elements based on composition
            num_total_atoms = len(atoms)
            symbols = []
            for element, fraction in self.config.composition.items():
                count = int(round(fraction * num_total_atoms))
                symbols.extend([element] * count)

            # Ensure the number of symbols matches the number of atoms
            while len(symbols) < num_total_atoms:
                symbols.append(random.choice(self.config.elements))
            symbols = symbols[:num_total_atoms]

            # Shuffle and assign to the atoms object
            random.shuffle(symbols)
            atoms.set_chemical_symbols(symbols)

            # Simple check for physical validity (avoids major overlaps)
            # A more robust check might be needed in a full implementation
            distances = atoms.get_all_distances(mic=True)
            min_dist = np.min(distances[np.nonzero(distances)])
            if min_dist < 1.0:
                # This is a fallback, in a real scenario we might regenerate or relax
                atoms.rattle(stdev=0.1, seed=42)

            structures.append(atoms)

        return structures
