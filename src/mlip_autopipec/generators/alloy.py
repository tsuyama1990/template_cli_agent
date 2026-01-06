from __future__ import annotations

import logging
import random
from typing import TYPE_CHECKING

import numpy as np
from ase.build import bulk

from .base import BaseStructureGenerator

if TYPE_CHECKING:
    from ase import Atoms
    from mlip_autopipec.config.models import SystemConfig

logger = logging.getLogger(__name__)

# Constants for alloy generation
DEFAULT_LATTICE_PARAM = 4.0  # Angstrom
MIN_SUPERCELL_ATOMS = 64


class AlloyGenerator(BaseStructureGenerator):
    """A generator for creating alloy structures."""

    def __init__(self, config: SystemConfig) -> None:
        """
        Initialise the AlloyGenerator.

        Args:
            config: The system configuration.
        """
        self.config = config

    def generate(self) -> list[Atoms]:
        """
        Generate a list of random solid-solution alloy structures.

        Returns:
            A list of ASE Atoms objects.
        """
        structures = []
        for _ in range(self.config.num_structures):
            primitive_cell = bulk(
                self.config.elements[0], self.config.lattice, a=DEFAULT_LATTICE_PARAM, cubic=True
            )
            supercell_size = self._get_supercell_size(primitive_cell)
            supercell = primitive_cell.repeat(supercell_size)  # type: ignore[no-untyped-call]
            self._set_composition(supercell)
            structures.append(supercell)

        logger.info(f"Generated {len(structures)} alloy structures.")
        return structures

    def _get_supercell_size(self, primitive_cell: Atoms) -> int:
        """Calculate the supercell size to get at least MIN_SUPERCELL_ATOMS."""
        num_atoms_primitive = len(primitive_cell)
        return int(np.ceil((MIN_SUPERCELL_ATOMS / num_atoms_primitive) ** (1 / 3)))

    def _set_composition(self, atoms: Atoms) -> None:
        """Set the composition of the supercell."""
        num_atoms = len(atoms)
        atomic_numbers = []
        for element, fraction in self.config.composition.items():
            count = round(num_atoms * fraction)
            atomic_numbers.extend([bulk(element).numbers[0]] * int(count))

        # Adjust for rounding errors
        while len(atomic_numbers) < num_atoms:
            atomic_numbers.append(bulk(self.config.elements[0]).numbers[0])
        atomic_numbers = atomic_numbers[:num_atoms]

        random.shuffle(atomic_numbers)
        atoms.set_atomic_numbers(atomic_numbers)  # type: ignore[no-untyped-call]
