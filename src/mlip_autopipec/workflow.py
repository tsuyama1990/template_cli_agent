import logging

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.orchestrator import Orchestrator

logger = logging.getLogger(__name__)


def run_cycle01_workflow(config: Cycle01Config):
    """
    Initializes dependencies and runs the full Cycle 01 workflow.

    This function sets up the necessary components (database wrapper, engines)
    and injects them into the Orchestrator, which then manages the execution
    of the labeling and training process.

    Args:
        config: The Cycle01Config object containing all settings.
    """
    try:
        # --- Dependency Injection Setup ---
        db_wrapper = AseDBWrapper(str(config.database_path))
        labeling_engine = LabelingEngine(config.dft_compute)
        training_engine = TrainingEngine(config.mlip_training)

        orchestrator = Orchestrator(
            db_wrapper=db_wrapper,
            labeling_engine=labeling_engine,
            training_engine=training_engine,
        )
        # --- End of Setup ---

        orchestrator.run_label_and_train_workflow()

    except Exception as e:
        logger.error(
            f"An unexpected error occurred during the workflow setup or execution: {e}",
            exc_info=True,
        )
        # Re-raise to allow the CLI to handle the exit
        raise
