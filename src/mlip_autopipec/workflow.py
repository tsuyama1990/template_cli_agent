from mlip_autopipec.config import DFTInputConfig, MLIPTrainingConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.modules.labeling_engine import LabelingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine


class WorkflowOrchestrator:
    """Orchestrates the entire MLIP generation workflow."""

    def __init__(
        self,
        dft_config: DFTInputConfig,
        training_config: MLIPTrainingConfig,
        db_path: str,
        qe_command: str,
    ):
        """Initializes the WorkflowOrchestrator.

        Args:
            dft_config: Configuration for the DFT calculations.
            training_config: Configuration for the MLIP training.
            db_path: Path to the ASE database.
            qe_command: The command to execute Quantum Espresso.
        """
        self.db_wrapper = AseDBWrapper(db_path)
        self.labeling_engine = LabelingEngine(
            dft_config, self.db_wrapper, qe_command
        )
        self.training_engine = TrainingEngine(training_config, self.db_wrapper)

    def run_labeling_for_id(self, id: int):
        """Runs the labeling process for a specific structure ID.

        Args:
            id: The ID of the structure to label.
        """
        print(f"Starting labeling for structure ID: {id}")
        self.labeling_engine.label_structure(id)
        print(f"Labeling complete for structure ID: {id}")

    def run_training(self):
        """Runs the MLIP model training process."""
        print("Starting training process.")
        try:
            self.training_engine.train()
            print("Training process finished.")
        except NotImplementedError as e:
            print(f"Training not implemented: {e}")
