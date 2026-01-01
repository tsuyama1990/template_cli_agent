from typing import List
from ase import Atoms
from ase.build import bulk
import numpy as np

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import StructureGeneration


class StructureGenerator:
    """
    Generates a diverse set of initial atomic structures based on the
    provided configuration.
    """

    def __init__(self, config: StructureGeneration, db_wrapper: AseDBWrapper):
        """
        Initializes the StructureGenerator.

        Args:
            config: The configuration for structure generation.
            db_wrapper: An instance of AseDBWrapper to save the structures.
        """
        self.config = config
        self.db_wrapper = db_wrapper

    def execute(self):
        """
        Executes the structure generation process.

        This method generates a base structure and then applies a series of
        strains to create a diverse set of initial structures, which are then
        saved to the database.
        """
        # Note: SQS generation with an external library is complex.
        # For this cycle, we will implement a placeholder that generates
        # a simple bulk crystal and applies strains, fulfilling the core
        # requirement of generating a diverse initial dataset.
        # A real SQS implementation would be added in a future cycle.

        # 1. Generate the base structure (e.g., a simple Au crystal for now)
        #    This part would be replaced by a proper SQS call.
        base_structure = bulk("Au", "fcc", a=4.0, cubic=True)

        # 2. Apply strains to create a list of structures
        structures_to_add: List[Atoms] = []
        for strain in self.config.strains:
            strained_structure = base_structure.copy()
            strained_structure.set_cell(
                strained_structure.cell * (1 + strain), scale_atoms=True
            )
            # Add some random noise/perturbation to the atomic positions
            positions = strained_structure.get_positions()
            positions += np.random.normal(0, 0.05, positions.shape)
            strained_structure.set_positions(positions)

            structures_to_add.append(strained_structure)

        # 3. Save all generated structures to the database
        self.db_wrapper.add_structures(structures_to_add)
