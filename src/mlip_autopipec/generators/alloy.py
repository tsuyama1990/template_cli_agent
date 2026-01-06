# src/mlip_autopipec/generators/alloy.py
"""
Implements the structure generator for multi-component alloys.
"""

import random

from ase import Atoms
from ase.build import bulk

from .base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """Generates random alloy structures."""

    def generate(self) -> list[Atoms]:
        """
        Creates a random alloy structure based on the system configuration.

        The method creates a supercell of the primary element and then
        randomly replaces atoms with other elements to match the desired
        composition.

        Returns:
            A list containing the single generated and validated alloy structure,
            or an empty list if validation fails.
        """
        system_cfg = self.config.system
        primary_element = system_cfg.elements[0]

        # Create a base crystal structure (fcc is a reasonable default)
        primitive_cell = bulk(primary_element, "fcc", a=4.0, cubic=True)

        # Create the supercell
        supercell = primitive_cell.repeat(system_cfg.supercell_size)  # type: ignore[no-untyped-call]
        n_atoms = len(supercell)

        # Determine the number of atoms of each element
        atom_counts = {el: round(comp * n_atoms) for el, comp in system_cfg.composition.items()}

        # Adjust counts to match total number of atoms due to rounding
        while sum(atom_counts.values()) < n_atoms:
            atom_counts[random.choice(list(atom_counts.keys()))] += 1  # noqa: S311
        while sum(atom_counts.values()) > n_atoms:
            atom_counts[random.choice(list(atom_counts.keys()))] -= 1  # noqa: S311

        # Create a list of atomic symbols with the desired composition
        new_symbols = []
        for element, count in atom_counts.items():
            new_symbols.extend([element] * count)

        random.shuffle(new_symbols)

        # Assign the new symbols to the supercell
        supercell.set_chemical_symbols(new_symbols)

        # Perform a small random rattle to break perfect symmetry
        supercell.rattle(stdev=0.1)

        if self._validate_structure(supercell):
            return [supercell]
        return []
