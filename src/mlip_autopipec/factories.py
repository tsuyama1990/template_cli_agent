from mlip_autopipec.config import (
    Settings,
)
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import IWorkflowOrchestrator
from mlip_autopipec.modules.labeling_engine import LabelingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.runners import SubprocessRunner
from mlip_autopipec.workflow import WorkflowOrchestrator


def create_workflow_orchestrator(
    db_path: str, settings: Settings
) -> IWorkflowOrchestrator:
    """
    Factory function to create a configured instance of the WorkflowOrchestrator.

    This function centralizes the creation of all dependencies, making the
    CLI and other entry points cleaner and more decoupled.

    Args:
        db_path: Path to the ASE database file.
        settings: The application settings.

    Returns:
        A fully configured object that implements the IWorkflowOrchestrator interface.
    """
    db_wrapper = AseDBWrapper(db_path)
    process_runner = SubprocessRunner()
    labeling_engine = LabelingEngine(
        settings.dft_input_configuration,
        process_runner,
        settings.qe_command,
    )
    training_engine = TrainingEngine(settings.mlip_training_configuration)

    return WorkflowOrchestrator(labeling_engine, training_engine, db_wrapper)
