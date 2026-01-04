import random

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell

from mlip_autopipec.common.pydantic_models import GenerationConfig, SystemConfig

from .base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    def __init__(self, system_config: SystemConfig, generation_config: GenerationConfig):
        self.system_config = system_config
        self.generation_config = generation_config

    def generate(self) -> list[Atoms]:
        structures = []
        for _ in range(self.system_config.num_initial_structures):
            # Create a base structure with the most common element
            base_element = max(self.system_config.composition, key=self.system_config.composition.get)
            primitive = bulk(base_element, "fcc", a=4.0, cubic=True)
            supercell = make_supercell(primitive, np.diag(self.system_config.supercell_size))

            # Randomly assign elements based on composition
            num_atoms = len(supercell)
            symbols = random.choices(
                list(self.system_config.composition.keys()),
                weights=list(self.system_config.composition.values()),
                k=num_atoms,
            )
            supercell.set_chemical_symbols(symbols)

            # Apply rattle and strain
            supercell.rattle(stdev=self.generation_config.rattle_std_dev)
            supercell.set_cell(
                supercell.get_cell() * (1 + self.generation_config.volumetric_strain),
                scale_atoms=True,
            )

            if self._validate_structure(supercell):
                structures.append(supercell)

        return structures

    def _validate_structure(self, atoms: Atoms) -> bool:
        """Checks for overlapping atoms."""
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        return np.min(distances) > self.generation_config.min_atomic_distance
