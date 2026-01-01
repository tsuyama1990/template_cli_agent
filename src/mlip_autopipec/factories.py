"""
This module contains factory functions for creating instances of the application's
core components.

Using factories helps to decouple the application's construction from its
business logic, adhering to the principles of Dependency Injection and
Inversion of Control. This makes the application more modular, easier to test,
and simpler to reconfigure.
"""

from mlip_autopipec.config import FullConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import IWorkflowOrchestrator
from mlip_autopipec.modules.labeling_engine import LabelingEngine
from mlip_autopipec.modules.structure_generator import StructureGenerator
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.runners import SubprocessRunner
from mlip_autopipec.workflow import WorkflowOrchestrator


def create_workflow_orchestrator(config: FullConfig) -> IWorkflowOrchestrator:
    """
    Factory function to create a fully configured `WorkflowOrchestrator`.

    This function is the central point for assembling the various components
    (like the database, labeling engine, and training engine) into a cohesive
    workflow orchestrator. It ensures that all dependencies are created and
    injected correctly based on the provided configuration.

    Args:
        config: The fully expanded and validated application configuration.

    Returns:
        A fully configured object that implements the `IWorkflowOrchestrator`
        interface.
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
    structure_generator = StructureGenerator(config, db_wrapper)

    return WorkflowOrchestrator(
        structure_generator, labeling_engine, training_engine, db_wrapper
    )
