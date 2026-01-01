from mlip_autopipec.config import FullConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import IWorkflowOrchestrator
from mlip_autopipec.modules.labeling_engine import LabelingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.runners import SubprocessRunner
from mlip_autopipec.workflow import WorkflowOrchestrator


def create_workflow_orchestrator(config: FullConfig) -> IWorkflowOrchestrator:
    """
    Factory function to create a configured instance of the WorkflowOrchestrator.

    This function centralizes the creation of all dependencies, making the
    CLI and other entry points cleaner and more decoupled.

    Args:
        config: The fully expanded application configuration.

    Returns:
        A fully configured object that implements the IWorkflowOrchestrator interface.
    """
    db_wrapper = AseDBWrapper(config.db_path)
    process_runner = SubprocessRunner()
    labeling_engine = LabelingEngine(
        config.dft_compute,  # Pass the DFTComputeConfig sub-model
        process_runner,
        config.qe_command,
    )
    training_engine = TrainingEngine(
        config.mlip_training  # Pass the MLIPTrainingConfig sub-model
    )

    return WorkflowOrchestrator(labeling_engine, training_engine, db_wrapper)
