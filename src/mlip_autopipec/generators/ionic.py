from __future__ import annotations

import math
import random
from typing import TYPE_CHECKING, Any

import numpy as np
from ase.build import bulk
from ase.data import oxidation_states

from mlip_autopipec.generators.base import BaseStructureGenerator

if TYPE_CHECKING:
    from ase import Atoms


class IonicGenerator(BaseStructureGenerator):
    """
    Generates charge-neutral ionic crystal structures.

    This generator determines the stoichiometry of a binary ionic compound
    based on the most common oxidation states of its constituent elements.
    It then constructs the material using a rocksalt (NaCl) crystal structure.
    """

    def __init__(self, config: dict[str, Any]) -> None:
        super().__init__(config)
        if len(self.elements) != 2:
            raise ValueError("IonicGenerator currently supports binary compounds only.")
        self.lattice_constant: float = config.get("lattice_constant", 4.0)

    def _get_stoichiometry(self) -> tuple[str, int, str, int]:
        """
        Calculate the charge-neutral stoichiometry for two elements.

        Returns:
            A tuple containing (cation_symbol, cation_count, anion_symbol, anion_count).
        """
        cation_symbol, anion_symbol = self.elements[0], self.elements[1]

        # Find the most common positive and negative oxidation states
        cation_ox_states = [s for s in oxidation_states[cation_symbol] if s > 0]
        anion_ox_states = [s for s in oxidation_states[anion_symbol] if s < 0]

        if not cation_ox_states or not anion_ox_states:
            raise ValueError(f"Could not determine oxidation states for {self.elements}")

        cation_charge = cation_ox_states[0]
        anion_charge = anion_ox_states[0]

        # Find the simplest integer ratio for charge neutrality
        common_divisor = math.gcd(cation_charge, abs(anion_charge))
        cation_count = abs(anion_charge) // common_divisor
        anion_count = cation_charge // common_divisor

        return cation_symbol, cation_count, anion_symbol, anion_count

    def generate(self) -> list[Atoms]:
        """
        Generate a list of ionic structures.

        Returns:
            A list of `ase.Atoms` objects with charge-neutral stoichiometry.
        """
        cat_sym, cat_n, an_sym, an_n = self._get_stoichiometry()
        symbols = [cat_sym] * cat_n + [an_sym] * an_n

        structures = []
        for _ in range(self.num_structures):
            # Create a disordered FCC structure as a starting point.
            # A more sophisticated generator would use specific crystal prototypes.
            num_atoms = len(symbols)

            # Create a supercell large enough to hold the atoms
            size = int(np.ceil(num_atoms**(1/3)))
            final_atoms = bulk(self.elements[0], "fcc", a=self.lattice_constant) * (size, size, size)

            # Ensure we have enough atoms in the supercell
            if len(final_atoms) < num_atoms:
                final_atoms = final_atoms * (2,1,1) # Double the size if needed

            # Truncate to the exact number of atoms required
            final_atoms = final_atoms[:num_atoms]

            # Set the chemical symbols according to the stoichiometry
            final_atoms.set_chemical_symbols(symbols)
            random.shuffle(final_atoms.symbols)

            final_atoms.rattle(stdev=0.1)
            structures.append(final_atoms)

        return structures
