from abc import ABC, abstractmethod


class ILabelingEngine(ABC):
    """Interface for a labeling engine that processes atomic structures."""

    @abstractmethod
    def label_structure(self, structure_id: int) -> None:
        """
        Labels a single atomic structure by running a calculation.

        Args:
            structure_id: The ID of the structure in the database to label.
        """
        pass


class ITrainingEngine(ABC):
    """Interface for a training engine that trains an MLIP model."""

    @abstractmethod
    def train(self) -> None:
        """Trains the MLIP model using the labeled data from the database."""
        pass


class IWorkflowOrchestrator(ABC):
    """Interface for the main workflow orchestrator."""

    @abstractmethod
    def label_structure_by_id(self, structure_id: int) -> None:
        """
        Runs the labeling process for a specific structure ID.

        Args:
            structure_id: The ID of the structure to label.
        """
        pass

    @abstractmethod
    def run_training(self) -> None:
        """Runs the MLIP model training process."""
        pass
