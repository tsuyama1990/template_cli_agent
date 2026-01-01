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
from mlip_autopipec.interfaces import IStructureGenerator, IWorkflowOrchestrator
from mlip_autopipec.modules.alloy_generator import AlloyGenerator
from mlip_autopipec.modules.labeling_engine import LabelingEngine
from mlip_autopipec.modules.molecule_generator import MoleculeGenerator
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.runners import SubprocessRunner
from mlip_autopipec.workflow import WorkflowOrchestrator


def create_structure_generator(config: FullConfig) -> IStructureGenerator:
    """
    Factory function to create the appropriate structure generator.

    This function acts as a dispatcher, reading the `structure_type` from the
    configuration and returning the correct generator implementation.

    Args:
        config: The fully expanded application configuration.

    Returns:
        An object that implements the `IStructureGenerator` interface.
    """
    structure_type = config.system.structure_type
    if structure_type == "alloy" or structure_type == "covalent":
        return AlloyGenerator(config)
    elif structure_type == "molecule":
        return MoleculeGenerator(config)
    else:
        # As more generators are added (e.g., for ionic), they will be added here.
        raise ValueError(f"Unsupported structure type: {structure_type}")


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
        config.dft_compute,
        process_runner,
        config.qe_command,
    )
    training_engine = TrainingEngine(config.mlip_training)
    structure_generator = create_structure_generator(config)

    return WorkflowOrchestrator(
        structure_generator, labeling_engine, training_engine, db_wrapper
    )
