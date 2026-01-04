import numpy as np
from ase import Atoms
from ase.build import bulk

from mlip_autopipec.common.pydantic_models import SystemConfig

from .base import BaseStructureGenerator


class AlloyGenerator(BaseStructureGenerator):
    """
    Concrete implementation for generating alloy structures.
    """

    def __init__(self, config: SystemConfig) -> None:
        self.config = config

    def generate(self, num_structures: int = 1) -> list[Atoms]:
        """
        Generates a list of random alloy structures.
        """
        structures = []
        for _ in range(num_structures):
            # Create a base structure using the first element
            primitive_cell = bulk(self.config.elements[0], "fcc", a=4.0, cubic=True)

            # Create a supercell
            supercell = primitive_cell * self.config.supercell_size

            # Get the number of atoms
            num_atoms = len(supercell)

            # Get the desired composition
            composition = self.config.composition

            # Create a list of atomic numbers based on composition
            elements = self.config.elements
            atomic_numbers = []
            for element, fraction in composition.items():
                count = round(num_atoms * fraction)
                atomic_numbers.extend([elements.index(element)] * count)

            # Ensure the list has the correct number of atoms
            while len(atomic_numbers) < num_atoms:
                atomic_numbers.append(elements.index(self.config.elements[0]))
            atomic_numbers = atomic_numbers[:num_atoms]

            # Shuffle the atomic numbers and assign to the supercell
            np.random.shuffle(atomic_numbers)
            supercell.set_atomic_numbers(atomic_numbers)

            # Apply a random rattle to the positions
            rattle_strength = 0.1
            supercell.rattle(stdev=rattle_strength, seed=np.random.randint(0, 100000))

            structures.append(supercell)

        return structures
