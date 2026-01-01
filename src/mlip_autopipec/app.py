from mlip_autopipec.config import Settings
from mlip_autopipec.factories import create_workflow_orchestrator
from mlip_autopipec.interfaces import IWorkflowOrchestrator


class Application:
    """The main application class."""

    def __init__(self, settings: Settings):
        """
        Initializes the application.

        Args:
            settings: The application settings.
        """
        self.settings = settings
        self.orchestrator: IWorkflowOrchestrator | None = None

    def _get_orchestrator(self, db_path: str | None) -> IWorkflowOrchestrator:
        """
        Gets a configured instance of the workflow orchestrator.

        Args:
            db_path: The path to the database file.

        Returns:
            A configured instance of the workflow orchestrator.
        """
        if self.orchestrator is None:
            db_path = db_path or self.settings.db_path
            self.orchestrator = create_workflow_orchestrator(db_path, self.settings)
        return self.orchestrator

    def label_structure(self, structure_id: int, db_path: str | None) -> None:
        """
        Labels a single atomic structure.

        Args:
            structure_id: The ID of the structure to label.
            db_path: The path to the database file.
        """
        orchestrator = self._get_orchestrator(db_path)
        orchestrator.label_structure_by_id(structure_id)

    def train(self, db_path: str | None) -> None:
        """
        Trains the MLIP model.

        Args:
            db_path: The path to the database file.
        """
        orchestrator = self._get_orchestrator(db_path)
        orchestrator.run_training()


def create_app() -> Application:
    """
    Creates a new application instance.

    Returns:
        A new application instance.
    """
    settings = Settings()
    return Application(settings)
