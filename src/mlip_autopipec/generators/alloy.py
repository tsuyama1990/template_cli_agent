import random

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.build.supercells import make_supercell

from mlip_autopipec.common.errors import PhysicsViolationError
from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.interfaces import BaseStructureGenerator

DEFAULT_LATTICE_CONSTANT = 3.5
MIN_INTERATOMIC_DISTANCE = 0.7  # In Angstrom

class AlloyGenerator(BaseStructureGenerator):
    """
    Generates initial alloy structures based on composition, ensuring physical validity.
    """

    def __init__(self, config: FullConfig):
        self.config = config.system
        self.num_structures = 1 # For now, generate one initial seed structure

    def generate(self) -> list[Atoms]:
        """
        Creates a list of physically valid, augmented alloy structures.

        Returns:
            A list containing a single generated ASE Atoms object.
        """
        atoms = self._create_random_alloy()

        # Apply augmentations
        atoms.rattle(stdev=0.1)

        # Apply volumetric strain (e.g., from -5% to +5%)
        strain = 1.0 + (random.random() - 0.5) * 0.1
        atoms.set_cell(atoms.cell * strain, scale_atoms=True)

        self._validate_structure(atoms)

        return [atoms]

    def _create_random_alloy(self) -> Atoms:
        """
        Builds a single random alloy structure.
        """
        primitive_cell = bulk(self.config.elements[0], 'sc', a=DEFAULT_LATTICE_CONSTANT)
        supercell_matrix = np.diag(self.config.supercell_size)
        supercell = make_supercell(primitive_cell, supercell_matrix)

        num_atoms = len(supercell)
        symbols = self._get_symbols_from_composition(num_atoms)
        random.shuffle(symbols)
        supercell.set_chemical_symbols(symbols)

        supercell.pbc = True
        return supercell

    def _get_symbols_from_composition(self, total_atoms: int) -> list[str]:
        """Calculates the number of atoms of each element based on composition."""
        symbols = []
        for element, fraction in self.config.composition.items():
            count = int(round(fraction * total_atoms))
            symbols.extend([element] * count)

        # Adjust for rounding errors to match total_atoms
        while len(symbols) < total_atoms:
            symbols.append(random.choice(list(self.config.composition.keys())))

        return symbols[:total_atoms]

    def _validate_structure(self, atoms: Atoms):
        """
        Checks if the structure is physically plausible.

        Raises:
            PhysicsViolationError: If atoms are too close together.
        """
        distances = atoms.get_all_distances(mic=True)
        # Ignore diagonal (distance to self)
        min_dist = np.min(distances[distances > 0])

        if min_dist < MIN_INTERATOMIC_DISTANCE:
            raise PhysicsViolationError(
                f"Generated structure has atoms too close: {min_dist:.2f} Å"
            )
