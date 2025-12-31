from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


class Orchestrator:
    """Orchestrates the MLIP-AutoPipe workflow."""

    def __init__(self, config: Cycle01Config):
        self.config = config
        self.db_wrapper = AseDBWrapper(str(config.database_path))

    def run_label_and_train_workflow(self):
        """Runs the Cycle 1 workflow: Label -> Train."""
        print("--- Starting Labeling Engine ---")
        labeling_engine = LabelingEngine(self.db_wrapper, self.config.dft_compute)
        labeling_engine.execute()
        print("--- Labeling Engine Finished ---")

        print("--- Starting Training Engine ---")
        training_engine = TrainingEngine(self.db_wrapper, self.config.mlip_training)
        training_engine.execute()
        print("--- Training Engine Finished ---")
