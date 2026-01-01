import logging

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine

logger = logging.getLogger(__name__)


class Orchestrator:
    """
    Orchestrates the MLIP-AutoPipe workflow for Cycle 1 by coordinating
    the labeling and training engines.
    """

    def __init__(
        self,
        db_wrapper: AseDBWrapper,
        labeling_engine: LabelingEngine,
        training_engine: TrainingEngine,
    ):
        """
        Initializes the Orchestrator with its dependencies.

        Args:
            db_wrapper: An instance of AseDBWrapper for database interactions.
            labeling_engine: An instance of LabelingEngine to run DFT calculations.
            training_engine: An instance of TrainingEngine to train the MLIP model.
        """
        self.db_wrapper = db_wrapper
        self.labeling_engine = labeling_engine
        self.training_engine = training_engine

    def run_label_and_train_workflow(self):
        """
        Executes the full Label -> Train workflow with error handling.
        """
        logger.info("--- Starting Cycle 01 Workflow: Label and Train ---")
        try:
            # 1. Run the Labeling Engine
            logger.info("Step 1: Running Labeling Engine...")
            rows_to_label = self.db_wrapper.get_rows_to_label()
            if not rows_to_label:
                logger.info("No new structures to label. Skipping labeling step.")
            else:
                structures_to_label = [
                    (row.id, row.toatoms()) for row in rows_to_label
                ]
                dft_results = self.labeling_engine.execute(structures_to_label)
                for row_id, results in dft_results:
                    self.db_wrapper.update_row_with_dft_results(row_id, results)
                logger.info("Labeling step completed successfully.")

            # 2. Run the Training Engine
            logger.info("Step 2: Running Training Engine...")
            labeled_rows = self.db_wrapper.get_all_labeled_rows()
            if not labeled_rows:
                logger.warning(
                    "No labeled structures found in the database. "
                    "Skipping training step."
                )
            else:
                atoms_to_train = [row.toatoms() for row in labeled_rows]
                self.training_engine.execute(atoms_to_train)
                logger.info("Training step completed successfully.")

            logger.info("--- Cycle 01 Workflow Finished Successfully ---")

        except Exception as e:
            logger.critical(
                f"Workflow failed due to a critical error: {e}", exc_info=True
            )
            # Re-raise the exception to be handled by the caller (e.g., the CLI)
            raise
