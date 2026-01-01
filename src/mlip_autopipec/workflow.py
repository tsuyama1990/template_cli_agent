# src/mlip_autopipec/workflow.py


from ase.io import read

from mlip_autopipec.configs.models import MainConfig
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


class WorkflowOrchestrator:
    """Orchestrates the entire MLIP generation workflow."""

    def __init__(self, config: MainConfig):
        """
        Initializes the WorkflowOrchestrator.

        Args:
            config: The main configuration object for the workflow.
        """
        self.config = config

    def run(self, database_path: str, input_file_path: str | None = None):
        """
        Executes the full labeling and training workflow.

        Args:
            database_path: The path to the ASE database file.
            input_file_path: Optional path to an XYZ file with initial structures.
        """
        print("Initializing Workflow")
        db_wrapper = AseDBWrapper(database_path)

        # Step 1: Populate database from input file if provided
        if input_file_path:
            print(f"Reading initial structures from {input_file_path}")
            try:
                initial_structures = read(input_file_path, index=":")
                if initial_structures:
                    db_wrapper.add_atoms(initial_structures)
                    print(
                        f"Added {len(initial_structures)} structures to the database."
                    )
            except FileNotFoundError:
                print(f"Error: Input file not found at {input_file_path}")
                return

        # Step 2: Run the Labeling Engine
        print("Starting Labeling Engine")
        labeling_engine = LabelingEngine(
            config=self.config.dft_compute, db_wrapper=db_wrapper
        )
        labeling_engine.run()
        print("Labeling Complete")

        # Step 3: Run the Training Engine
        print("Starting Training Engine")
        training_engine = TrainingEngine(
            config=self.config.mlip_training, db_wrapper=db_wrapper
        )
        training_engine.run()
        print("Training Complete")
