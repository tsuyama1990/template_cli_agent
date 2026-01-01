"""
This module defines the main `Application` class that orchestrates the workflow.

It acts as the primary entry point for the business logic, decoupling the CLI
from the core implementation of the labeling and training processes.
"""

from mlip_autopipec.config import FullConfig
from mlip_autopipec.factories import create_workflow_orchestrator
from mlip_autopipec.interfaces import IWorkflowOrchestrator


class Application:
    """
    The main application class that encapsulates the workflow logic.

    This class is initialized with a complete configuration object and provides
    methods to execute the main tasks of the application, such as labeling
    structures and training the MLIP model.
    """

    def __init__(self, config: FullConfig):
        """
        Initializes the application with a full configuration.

        Args:
            config: The application's full, expanded configuration object.
        """
        self.config = config
        self.orchestrator: IWorkflowOrchestrator | None = None

    def _get_orchestrator(self) -> IWorkflowOrchestrator:
        """
        Lazily initializes and returns a configured workflow orchestrator.

        This method ensures that the orchestrator and its dependencies (like the
        database connection) are only created when a workflow task is actually
        called.

        Returns:
            A fully configured instance of the workflow orchestrator.
        """
        if self.orchestrator is None:
            self.orchestrator = create_workflow_orchestrator(self.config)
        return self.orchestrator

    def label_structure(self, structure_id: int) -> None:
        """
        Labels a single atomic structure by its ID.

        Args:
            structure_id: The ID of the structure in the database to label.
        """
        orchestrator = self._get_orchestrator()
        orchestrator.label_structure_by_id(structure_id)

    def train(self) -> None:
        """
        Trains the MLIP model using all available labeled data.
        """
        orchestrator = self._get_orchestrator()
        orchestrator.run_training()

    def run(self) -> None:
        """
        Runs the full MLIP generation workflow.
        """
        orchestrator = self._get_orchestrator()
        orchestrator.run()
