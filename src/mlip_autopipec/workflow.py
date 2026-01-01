from mlip_autopipec.interfaces import (
    ILabelingEngine,
    ITrainingEngine,
    IWorkflowOrchestrator,
)


class WorkflowOrchestrator(IWorkflowOrchestrator):
    """
    Orchestrates the entire MLIP generation workflow.

    This class is responsible for coordinating the different components of the
    pipeline, such as the labeling and training engines, to execute the
    complete workflow for generating an MLIP.
    """

    def __init__(
        self,
        labeling_engine: ILabelingEngine,
        training_engine: ITrainingEngine,
    ):
        """
        Initializes the WorkflowOrchestrator.

        Args:
            labeling_engine: An object that implements the ILabelingEngine interface.
            training_engine: An object that implements the ITrainingEngine interface.
        """
        self.labeling_engine = labeling_engine
        self.training_engine = training_engine

    def label_structure_by_id(self, structure_id: int) -> None:
        """
        Runs the labeling process for a specific structure ID.

        Args:
            structure_id: The ID of the structure to label.
        """
        print(f"Starting labeling for structure ID: {structure_id}")
        self.labeling_engine.label_structure(structure_id)
        print(f"Labeling complete for structure ID: {structure_id}")

    def run_training(self) -> None:
        """Runs the MLIP model training process."""
        print("Starting training process.")
        try:
            self.training_engine.train()
            print("Training process finished.")
        except NotImplementedError as e:
            print(f"Training not implemented: {e}")
