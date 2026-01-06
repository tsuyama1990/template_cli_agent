# src/mlip_autopipec/generators/alloy.py
"""Concrete implementation of a structure generator for alloys."""

import random

from ase import Atoms
from ase.build import bulk

from mlip_autopipec.common.pydantic_models import SystemConfig
from mlip_autopipec.generators.base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """Generates random alloy structures based on the provided configuration."""

    def __init__(self, config: SystemConfig) -> None:
        self.config = config

    def generate(self) -> list[Atoms]:
        """
        Generates a single random alloy structure.

        The current implementation creates a primitive cell of the first element,
        builds a supercell, and then randomly replaces atoms to match the
        desired composition.

        Returns:
            A list containing a single ASE Atoms object for the generated alloy.
        """
        # Create a base structure using the first element in the list
        base_element = self.config.elements[0]
        atoms = bulk(base_element, "fcc", a=4.0, cubic=True)  # type: ignore[no-untyped-call]

        # Create the supercell
        atoms = atoms.repeat(self.config.supercell_size)  # type: ignore[no-untyped-call]

        # Get the total number of atoms
        n_atoms = len(atoms)

        # Create a list of atoms to be placed
        target_atoms = []
        for element, fraction in self.config.composition.items():
            count = round(n_atoms * fraction)
            target_atoms.extend([element] * count)

        # Fill any rounding errors with the base element
        while len(target_atoms) < n_atoms:
            target_atoms.append(base_element)

        # Truncate if we have too many
        target_atoms = target_atoms[:n_atoms]

        # Shuffle and assign to the atoms object
        random.shuffle(target_atoms)
        atoms.set_chemical_symbols(target_atoms)

        # Apply a random rattle to break symmetry
        atoms.rattle(stdev=0.1)

        return [atoms]
