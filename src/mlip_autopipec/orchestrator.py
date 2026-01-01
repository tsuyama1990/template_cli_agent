import logging

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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
        logger.info("--- Starting Cycle 01 Workflow: Label and Train ---")

        # 1. Initialize and run the Labeling Engine
        labeling_engine = LabelingEngine(self.config.dft_compute)
        rows_to_label = self.db_wrapper.get_rows_to_label()
        structures_to_label = [(row.id, row.toatoms()) for row in rows_to_label]
        dft_results = labeling_engine.execute(structures_to_label)
        for row_id, results in dft_results:
            self.db_wrapper.update_row_with_dft_results(row_id, results)

        # 2. Initialize and run the Training Engine
        training_engine = TrainingEngine(self.config.mlip_training)
        labeled_rows = self.db_wrapper.get_all_labeled_rows()
        atoms_to_train = [row.toatoms() for row in labeled_rows]
        training_engine.execute(atoms_to_train)

        logger.info("--- Cycle 01 Workflow Finished ---")
