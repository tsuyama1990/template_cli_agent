from mlip_autopipec.config import FullConfig
from mlip_autopipec.factories import create_workflow_orchestrator
from mlip_autopipec.interfaces import IWorkflowOrchestrator


class Application:
    """The main application class."""

    def __init__(self, config: FullConfig):
        """
        Initializes the application.

        Args:
            config: The application's full, expanded configuration.
        """
        self.config = config
        self.orchestrator: IWorkflowOrchestrator | None = None

    def _get_orchestrator(self) -> IWorkflowOrchestrator:
        """
        Gets a configured instance of the workflow orchestrator.

        Returns:
            A configured instance of the workflow orchestrator.
        """
        if self.orchestrator is None:
            self.orchestrator = create_workflow_orchestrator(self.config)
        return self.orchestrator

    def label_structure(self, structure_id: int) -> None:
        """
        Labels a single atomic structure.

        Args:
            structure_id: The ID of the structure to label.
        """
        orchestrator = self._get_orchestrator()
        orchestrator.label_structure_by_id(structure_id)

    def train(self) -> None:
        """
        Trains the MLIP model.
        """
        orchestrator = self._get_orchestrator()
        orchestrator.run_training()
