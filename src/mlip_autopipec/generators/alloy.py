import random

import numpy as np
from ase import Atoms
from ase.build import bulk

from mlip_autopipec.generators.base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """A generator for creating simple alloy structures."""

    def generate(self) -> list[Atoms]:
        """
        Generates a random alloy structure based on the configuration.

        Note: This is a simplified implementation for Cycle 1.
        """
        # Create a base crystal structure (e.g., BCC for Fe-based alloy)
        # A more robust implementation would select this based on elements.
        primitive_cell = bulk(self.config.elements[0], "bcc", a=2.87, cubic=True)

        # Create a supercell
        supercell = primitive_cell * tuple(self.config.supercell_size)

        # Get the total number of atoms
        num_atoms = len(supercell)

        # Assign elements based on composition
        symbols = []
        remaining_atoms = num_atoms
        for element, fraction in self.config.composition.items():
            count = round(fraction * num_atoms)
            symbols.extend([element] * count)
            remaining_atoms -= count

        # Add remaining atoms due to rounding to the most numerous element
        if remaining_atoms > 0:
            most_common_element = max(
                self.config.composition, key=lambda k: self.config.composition[k]
            )
            symbols.extend([most_common_element] * remaining_atoms)

        # Shuffle the symbols and assign them to the supercell
        random.shuffle(symbols)
        supercell.set_chemical_symbols(symbols)

        # Apply a random rattle to break symmetry
        positions = supercell.get_positions()
        noise = np.random.uniform(-0.1, 0.1, size=positions.shape)
        supercell.set_positions(positions + noise)
        supercell.wrap()

        return [supercell]
