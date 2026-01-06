"""A concrete implementation for generating alloy structures."""

import random
import numpy as np
from typing import List
from ase import Atoms
from ase.build import bulk
from ase.build.supercells import make_supercell
import logging

from mlip_autopipec.generators.base import BaseStructureGenerator
from mlip_autopipec.common.pydantic_models import SystemConfig

# Configure logger for this module
logger = logging.getLogger(__name__)

class AlloyGenerator(BaseStructureGenerator):
    """
    Generates initial structures for a multi-component alloy.

    This generator creates a supercell and randomly assigns atomic species
    to the lattice sites according to the specified composition.
    """

    def __init__(self, config: SystemConfig) -> None:
        """
        Initializes the AlloyGenerator.

        Args:
            config: The system configuration Pydantic model.
        """
        self.config = config
        logger.info(f"Initialized AlloyGenerator with config: {config}")

    def generate(self) -> List[Atoms]:
        """
        Generates a single random alloy structure.

        For Cycle 1, this generates one structure. It can be extended to
        generate multiple distinct structures if needed.

        Returns:
            A list containing a single ASE Atoms object for the alloy.
        """
        logger.info("Generating alloy structure...")

        # 1. Create a primitive cell of the primary element
        # Using a simple crystal structure as a base. 'fcc' is a reasonable default.
        primary_element = self.config.elements[0]
        primitive_cell = bulk(primary_element, 'fcc', a=4.0, cubic=True)

        # 2. Create the supercell transformation matrix
        supercell_matrix = np.diag(self.config.supercell_size)

        # 3. Build the supercell
        supercell = make_supercell(primitive_cell, supercell_matrix) # type: ignore
        n_atoms = len(supercell)
        logger.info(f"Created a supercell with {n_atoms} atoms.")

        # 4. Determine the number of atoms for each element
        symbols = []
        remaining_atoms = n_atoms
        for i, element in enumerate(self.config.elements):
            if i == len(self.config.elements) - 1:
                # Assign the rest to the last element to avoid rounding errors
                count = remaining_atoms
            else:
                count = round(n_atoms * self.config.composition[element])

            symbols.extend([element] * count)
            remaining_atoms -= count

        if len(symbols) != n_atoms:
            # This should ideally not happen with the rounding logic, but as a safeguard:
            msg = "Mismatch in atom counts during symbol assignment."
            raise ValueError(msg)

        # 5. Randomly assign symbols to atomic positions
        random.shuffle(symbols)
        supercell.set_chemical_symbols(symbols)

        # 6. Apply a small random rattle to break symmetry
        # This is good practice to start MD simulations.
        supercell.rattle(stdev=0.01, seed=42)
        logger.info("Applied random rattle to the structure.")

        # Perform a simple validation check
        if self._has_overlapping_atoms(supercell):
            logger.warning("Generated structure has overlapping atoms. This may cause issues.")

        logger.info("Successfully generated alloy structure.")
        return [supercell]

    def _has_overlapping_atoms(self, atoms: Atoms, threshold: float = 0.7) -> bool:
        """
        Checks if any two atoms in the structure are closer than a given threshold.

        Args:
            atoms: The ASE Atoms object to check.
            threshold: The minimum allowed distance between atoms in Angstroms.

        Returns:
            True if overlapping atoms are found, False otherwise.
        """
        distances = atoms.get_all_distances(mic=True) # type: ignore
        # Set diagonal to a large value to ignore self-distances
        np.fill_diagonal(distances, np.inf)
        min_distance = np.min(distances)

        if min_distance < threshold:
            logger.warning(
                f"Overlap detected: Minimum inter-atomic distance is {min_distance:.2f} Å, which is below the threshold of {threshold} Å."
            )
            return True
        return False
