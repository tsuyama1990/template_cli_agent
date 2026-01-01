from mlip_autopipec.config import (
    DFTInputConfig,
    MLIPTrainingConfig,
    ModelType,
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
    dft_config = DFTInputConfig(
        pseudopotentials={"H": "H.UPF", "O": "O.UPF"},
        kpoints=(1, 1, 1),
        ecutwfc=60,
        control={"calculation": "scf"},
    )

    training_config = MLIPTrainingConfig(
        model_type=ModelType.ACE,
        r_cut=5.0,
        loss_weights={"energy": 1.0, "forces": 100.0},
    )

    db_wrapper = AseDBWrapper(db_path)
    process_runner = SubprocessRunner()
    labeling_engine = LabelingEngine(
        dft_config, db_wrapper, process_runner, settings.qe_command
    )
    training_engine = TrainingEngine(training_config, db_wrapper)

    return WorkflowOrchestrator(labeling_engine, training_engine)
