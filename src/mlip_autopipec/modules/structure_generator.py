"""Module A: Generates initial atomic structures."""
import logging

from ase import Atoms

# Local imports for type hinting. Will cause circular import if used directly.
if hasattr(__import__("typing"), "TYPE_CHECKING") and __import__("typing").TYPE_CHECKING:
    from mlip_autopipec.config import FullConfig
    from mlip_autopipec.database import AseDBWrapper


logger = logging.getLogger(__name__)


class StructureGenerator:
    """
    Generates a diverse and physically relevant set of initial atomic structures
    based on the chemical composition and material type.
    """

    def __init__(self, config: "FullConfig", db_wrapper: "AseDBWrapper"):
        """
        Initializes the StructureGenerator.

        Args:
            config: The full, expanded configuration object.
            db_wrapper: An instance of the ASE database wrapper.
        """
        self.config = config
        self.db = db_wrapper
        self.system_info = config.system

    def generate(self):
        """
        Main entry point for structure generation.

        Checks if the database is empty. If so, it proceeds to generate
        and save the initial set of structures.
        """
        if self.db.is_empty():
            logger.info("Database is empty. Generating initial structures...")
            structures = self._dispatcher()
            for atoms in structures:
                self.db.add_atoms(atoms, state="unlabeled")
            message = f"Successfully generated and saved {len(structures)} new structures."
            logger.info(message)
            print(message)  # Also print to stdout for CLI visibility
        else:
            logger.info("Database is not empty. Skipping initial structure generation.")

    def _dispatcher(self) -> list[Atoms]:
        """
        Dispatches to the correct generation method based on structure_type.
        """
        structure_type = self.system_info.structure_type
        logger.info(f"Structure type '{structure_type}' detected. Dispatching to appropriate generator.")

        if structure_type == "alloy":
            return self._generate_sqs_for_alloy()
        elif structure_type == "molecule":
            return self._generate_nms_for_molecule()
        elif structure_type == "ionic":
            return self._generate_airss_for_ionic()
        elif structure_type == "covalent":
            return self._generate_rattle_for_covalent()
        else:
            raise ValueError(f"Unknown structure type: {structure_type}")

    def _generate_sqs_for_alloy(self) -> list[Atoms]:
        """Generates Special Quasirandom Structures (SQS) for alloys."""
        import numpy as np

        # Simple heuristic to create a base structure
        # Alternative implementation using ASE to avoid external binary dependencies.
        # This creates a randomized supercell, which is a good approximation for an initial guess.
        from ase.build import bulk
        from ase.build.supercells import make_supercell
        from pymatgen.core import Composition

        # Create a base primitive cell (e.g., simple cubic)
        # We use a placeholder element that will be replaced.
        primitive_cell = bulk(self.system_info.elements[0], "sc", a=3.5, cubic=True)

        # Define scaling to get a supercell of reasonable size
        target_atoms = 64
        scaling_factor = int(np.ceil((target_atoms / len(primitive_cell)) ** (1 / 3)))
        scaling_matrix = np.diag([scaling_factor] * 3)
        supercell = make_supercell(primitive_cell, scaling_matrix)

        # Determine the number of atoms of each type based on composition
        composition = Composition(self.system_info.composition).fractional_composition
        total_atoms = len(supercell)
        atomic_numbers = []
        for element, fraction in composition.items():
            count = round(fraction * total_atoms)
            atomic_numbers.extend([element.Z] * int(count))

        # Ensure the total count matches the supercell size
        while len(atomic_numbers) < total_atoms:
            atomic_numbers.append(composition.elements[-1].Z)
        atomic_numbers = atomic_numbers[:total_atoms]

        # Randomly assign the elements to the sites
        np.random.shuffle(atomic_numbers)
        supercell.set_atomic_numbers(atomic_numbers)

        base_atoms = supercell
        generated_structures = [base_atoms]

        # Generate strained structures
        for _ in range(self.config.simulation.initial_structures_to_generate - 1):
            strained_atoms = base_atoms.copy()
            strain = np.random.uniform(-0.05, 0.05, (3, 3))
            strain = (strain + strain.T) / 2  # Symmetrize the strain tensor
            new_cell = strained_atoms.get_cell() @ (np.identity(3) + strain)
            strained_atoms.set_cell(new_cell, scale_atoms=True)
            generated_structures.append(strained_atoms)

        return generated_structures

    def _generate_nms_for_molecule(self) -> list[Atoms]:
        """Generates distorted structures using Normal Mode Sampling (NMS)."""
        import warnings
        warnings.warn("NMS generation is not yet implemented", UserWarning)
        # TODO: Implement NMS logic.
        return []

    def _generate_airss_for_ionic(self) -> list[Atoms]:
        """Generates structures using Ab Initio Random Structure Searching (AIRSS)."""
        logger.warning("AIRSS generation is not yet implemented. Returning empty list.")
        # TODO: Implement AIRSS logic.
        return []

    def _generate_rattle_for_covalent(self) -> list[Atoms]:
        """Generates distorted structures by 'rattling' atoms."""
        logger.warning("Rattle generation is not yet implemented. Returning empty list.")
        # TODO: Implement rattling logic.
        return []
