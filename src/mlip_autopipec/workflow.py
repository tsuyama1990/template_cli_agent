import logging

from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.orchestrator import Orchestrator

logger = logging.getLogger(__name__)


def run_cycle01_workflow(config: Cycle01Config):
    """
    Runs the full Cycle 01 workflow.

    Args:
        config: The Cycle01Config object.
    """
    try:
        orchestrator = Orchestrator(config)
        orchestrator.run_label_and_train_workflow()
    except Exception as e:
        logger.error(f"An unexpected error occurred during the workflow: {e}")
        raise
