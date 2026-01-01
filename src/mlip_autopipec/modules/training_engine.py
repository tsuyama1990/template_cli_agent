from ase import Atoms

from mlip_autopipec.config import DFTResult, MLIPTrainingConfig
from mlip_autopipec.interfaces import ITrainingEngine


class TrainingEngine(ITrainingEngine):
    """
    Handles the training of the MLIP model.

    This class takes labeled data, prepares it for the training library
    (e.g., pacemaker), and executes the training process to produce an
    MLIP model artifact.
    """

    def __init__(
        self,
        mlip_training_configuration: MLIPTrainingConfig,
    ):
        """
        Initializes the TrainingEngine.

        Args:
            mlip_training_configuration: Configuration for the MLIP training.
        """
        self.mlip_training_configuration = mlip_training_configuration

    def train(self, training_data: list[tuple[Atoms, DFTResult]]) -> None:
        """
        Trains the MLIP model using the provided labeled data.

        Args:
            training_data: A list of tuples, each containing an Atoms object
                           and its corresponding DFTResult.
        """
        if not training_data:
            print("No labeled data provided to train on.")
            return

        print(f"Received {len(training_data)} labeled structures for training.")
        print(f"Training model of type: {self.mlip_training_configuration.model_type.value}")

        raise NotImplementedError("Training with pacemaker-python is not yet implemented.")
