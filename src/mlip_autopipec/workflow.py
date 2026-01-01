from .database import AseDBWrapper
from .modules.labeling_engine import LabelingEngine
from .modules.training_engine import TrainingEngine
from .config import DFTInputConfig, MLIPTrainingConfig
from .interfaces import IWorkflowOrchestrator
import yaml

class WorkflowOrchestrator(IWorkflowOrchestrator):
    """Orchestrates the entire MLIP generation workflow."""

    def __init__(self, config_path: str):
        """
        Initializes the orchestrator with a configuration file.

        Args:
            config_path: Path to the main YAML configuration file.
        """
        self.config = self._load_config(config_path)

        # In a larger app, this might be handled by a DI framework
        db_path = self.config.get("database_path", "asedb.db")
        self.db_wrapper = AseDBWrapper(db_path)

        dft_config = DFTInputConfig(**self.config["dft_compute"])
        qe_command = self.config["qe_command"]
        self.labeling_engine = LabelingEngine(self.db_wrapper, dft_config, qe_command)

        training_config = MLIPTrainingConfig(**self.config["mlip_training"])
        self.training_engine = TrainingEngine(self.db_wrapper, training_config)

    def _load_config(self, config_path: str) -> dict:
        """Loads the YAML configuration file."""
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)

    def run_labeling_for_id(self, uid: int):
        """
        Runs the labeling process for a single structure.

        Args:
            uid: The ID of the structure to label.
        """
        self.labeling_engine.label_structure(uid)

    def run_training(self):
        """Runs the MLIP model training process."""
        self.training_engine.train()

def create_workflow_orchestrator(config_path: str) -> IWorkflowOrchestrator:
    """
    Factory function to create an instance of the workflow orchestrator.

    This provides a single point of construction and decouples the client
    (e.g., the CLI) from the concrete implementation details.

    Args:
        config_path: Path to the main YAML configuration file.

    Returns:
        An object that conforms to the IWorkflowOrchestrator interface.
    """
    return WorkflowOrchestrator(config_path)
