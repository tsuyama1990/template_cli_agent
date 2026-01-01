from abc import ABC, abstractmethod

class IWorkflowOrchestrator(ABC):
    """
    Defines the public interface for the workflow orchestrator.
    This abstraction decouples the CLI from the concrete implementation.
    """

    @abstractmethod
    def run_labeling_for_id(self, uid: int):
        """
        Runs the labeling process for a single structure.

        Args:
            uid: The ID of the structure to label.
        """
        pass

    @abstractmethod
    def run_training(self):
        """Runs the MLIP model training process."""
        pass
