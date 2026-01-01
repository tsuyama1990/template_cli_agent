import yaml
import ase.io

from mlip_autopipec.configs.models import MainConfig
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


class WorkflowOrchestrator:
    """
    Orchestrates the entire MLIP generation workflow, from data ingestion to training.
    """

    def __init__(
        self,
        config_path: str,
        db_path: str,
        input_file_path: str,
    ):
        print("Initializing Workflow...")
        self.config = self._load_config(config_path)
        self.db_path = db_path
        self.input_file_path = input_file_path
        self.db_wrapper = AseDBWrapper()

    def _load_config(self, path: str) -> MainConfig:
        """Loads and validates the main configuration file."""
        print(f"Loading configuration from {path}")
        with open(path, "r") as f:
            config_dict = yaml.safe_load(f)
        return MainConfig(**config_dict)

    def _ingest_initial_structures(self):
        """Reads initial structures from an XYZ file and adds them to the database."""
        if self.input_file_path:
            print(f"Ingesting initial structures from {self.input_file_path}")
            initial_structures = ase.io.read(self.input_file_path, index=":")
            self.db_wrapper.add_atoms(initial_structures)
            print(f"Added {len(initial_structures)} structures to the database.")

    def run_workflow(self):
        """
        Executes the full labeling and training workflow in sequence.
        """
        # 1. Connect to the database
        self.db_wrapper.connect(self.db_path)
        print(f"Connected to database at {self.db_path}")

        # 2. Ingest initial data if provided
        self._ingest_initial_structures()

        # 3. Instantiate and run the Labeling Engine
        labeling_engine = LabelingEngine(
            config=self.config.dft_compute,
            db_wrapper=self.db_wrapper,
        )
        labeling_engine.run()

        # 4. Instantiate and run the Training Engine
        training_engine = TrainingEngine(
            config=self.config.mlip_training,
            db_wrapper=self.db_wrapper,
        )
        training_engine.run()

        print("Workflow finished successfully.")
