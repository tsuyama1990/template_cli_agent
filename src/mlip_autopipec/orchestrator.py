from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine

class Orchestrator:
    """
    Orchestrates the MLIP-AutoPipe workflow for Cycle 1.
    """

    def __init__(self, config: Cycle01Config):
        """
        Initializes the Orchestrator.

        Args:
            config: The Cycle01Config object containing all settings.
        """
        self.config = config
        self.db_wrapper = AseDBWrapper(str(config.database_path))

    def run_label_and_train_workflow(self):
        """
        Executes the full Label -> Train workflow.
        """
        print("--- Starting Cycle 01 Workflow: Label and Train ---")

        # 1. Initialize and run the Labeling Engine
        labeling_engine = LabelingEngine(self.config.dft_compute, self.db_wrapper)
        labeling_engine.execute()

        # 2. Initialize and run the Training Engine
        training_engine = TrainingEngine(self.config.mlip_training, self.db_wrapper)
        training_engine.execute()

        print("--- Cycle 01 Workflow Finished ---")
