from mlip_autopipec.config import MLIPTrainingConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import ITrainingEngine


class TrainingEngine(ITrainingEngine):
    """
    Handles the training of the MLIP model.

    This class queries the database for labeled data, prepares it for the
    training library (e.g., pacemaker), and executes the training process
    to produce an MLIP model artifact.
    """

    def __init__(
        self,
        training_config: MLIPTrainingConfig,
        db_wrapper: AseDBWrapper,
    ):
        """
        Initializes the TrainingEngine.

        Args:
            training_config: Configuration for the MLIP training.
            db_wrapper: Wrapper for the ASE database.
        """
        self.training_config = training_config
        self.db_wrapper = db_wrapper

    def train(self) -> None:
        """Trains the MLIP model using the labeled data from the database."""
        labeled_data = self.db_wrapper.get_all_labeled_atoms()
        if not labeled_data:
            print("No labeled data found to train on.")
            return

        # The actual training logic using pacemaker will be implemented in a future cycle.
        # For now, we'll just print a message.
        print(f"Found {len(labeled_data)} labeled structures for training.")
        print(f"Training model of type: {self.training_config.model_type.value}")

        raise NotImplementedError(
            "Training with pacemaker-python is not yet implemented."
        )
