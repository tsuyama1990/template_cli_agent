

from .config.models import FullConfig
from .data.database import AseDB
from .modules.c_labelling_engine import LabellingEngine
from .modules.d_training_engine import TrainingEngine


class Orchestrator:
    """The central coordinator for the MLIP-AutoPipe workflow."""

    def __init__(self, config: FullConfig):
        """
        Initializes the Orchestrator with a validated configuration.

        Args:
            config: The full, validated configuration object.
        """
        self.config = config
        self.db = AseDB(db_path=config.ase_db_path)
        self.labelling_engine = LabellingEngine(dft_config=config.dft_compute)
        self.training_engine = TrainingEngine(train_config=config.training)

    def run_cycle_01(self):
        """
        Executes the full workflow for Cycle 01:
        1. Finds structures that need labelling.
        2. Runs the Labelling Engine on them.
        3. Updates the database with the results.
        4. Gathers all labelled data.
        5. Runs the Training Engine.
        """
        print("--- Starting Cycle 01 Workflow ---")

        # 1. Find and label structures
        structures_to_label = self.db.get_structures_by_status("needs_labelling")
        if not structures_to_label:
            print("No structures found with status 'needs_labelling'.")
        else:
            print(f"Found {len(structures_to_label)} structures to label.")
            for i, atoms in enumerate(structures_to_label):
                print(f"  Labelling structure {i+1}/{len(structures_to_label)}...")
                try:
                    labelled_atoms = self.labelling_engine.run(atoms)
                    self.db.update_structure(labelled_atoms, status="labelled")
                    print(f"  Structure {i+1} labelled successfully.")
                except Exception as e:
                    print(f"  Error labelling structure {i+1}: {e}")
                    # In a real scenario, we might want more robust ID tracking
                    self.db.update_structure(atoms, status="failed_labelling")

        # 2. Gather data and train model
        labelled_structures = self.db.get_structures_by_status("labelled")
        if not labelled_structures:
            print("No labelled structures found. Skipping training.")
            return

        print(f"\nFound {len(labelled_structures)} labelled structures for training.")
        self.training_engine.run(labelled_structures)

        print("\n--- Cycle 01 Workflow Finished ---")
