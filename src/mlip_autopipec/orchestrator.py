import logging

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.modules.a_structure_generator import StructureGenerator
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine

logger = logging.getLogger(__name__)


class Orchestrator:
    """
    Orchestrates the MLIP-AutoPipe workflow by coordinating the various
    engine modules.
    """

    def __init__(
        self,
        db_wrapper: AseDBWrapper,
        structure_generator: StructureGenerator,
        labeling_engine: LabelingEngine,
        training_engine: TrainingEngine,
    ):
        """
        Initializes the Orchestrator with its dependencies.
        """
        self.db_wrapper = db_wrapper
        self.structure_generator = structure_generator
        self.labeling_engine = labeling_engine
        self.training_engine = training_engine

    def run_full_pipeline(self):
        """
        Executes the full Generate -> Label -> Train workflow with error handling.
        """
        logger.info("--- Starting Full Workflow: Generate, Label, and Train ---")
        try:
            # 1. Run the Structure Generator
            logger.info("Step 1: Running Structure Generator...")
            self.structure_generator.execute()
            logger.info("Structure generation completed successfully.")

            # 2. Run the Labeling Engine
            logger.info("Step 2: Running Labeling Engine...")
            rows_to_label = self.db_wrapper.get_rows_to_label()
            if not rows_to_label:
                logger.warning("No new structures to label. Skipping labeling step.")
            else:
                structures_to_label = [
                    (row.id, row.toatoms()) for row in rows_to_label
                ]
                dft_results = self.labeling_engine.execute(structures_to_label)
                for row_id, results in dft_results:
                    self.db_wrapper.update_row_with_dft_results(row_id, results)
                logger.info("Labeling step completed successfully.")

            # 3. Run the Training Engine
            logger.info("Step 3: Running Training Engine...")
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

            logger.info("--- Full Workflow Finished Successfully ---")

        except Exception as e:
            logger.critical(
                f"Workflow failed due to a critical error: {e}", exc_info=True
            )
            # Re-raise the exception to be handled by the caller (e.g., the CLI)
            raise
