from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import (
    ILabelingEngine,
    IStructureGenerator,
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
        structure_generator: IStructureGenerator,
        labeling_engine: ILabelingEngine,
        training_engine: ITrainingEngine,
        db_wrapper: AseDBWrapper,
    ):
        """
        Initializes the WorkflowOrchestrator.

        Args:
            structure_generator: An object that implements the IStructureGenerator interface.
            labeling_engine: An object that implements the ILabelingEngine interface.
            training_engine: An object that implements the ITrainingEngine interface.
            db_wrapper: A wrapper for the ASE database.
        """
        self.structure_generator = structure_generator
        self.labeling_engine = labeling_engine
        self.training_engine = training_engine
        self.db_wrapper = db_wrapper

    def run(self) -> None:
        """Runs the main workflow."""
        if self.db_wrapper.is_empty():
            print("Database is empty. Generating initial structures...")
            structures = self.structure_generator.generate()
            for atoms in structures:
                self.db_wrapper.add_atoms(atoms, state="unlabeled")
            print(f"Successfully saved {len(structures)} new structures to the database.")
        else:
            print("Database is not empty. Skipping initial structure generation.")

        # The rest of the workflow logic will be added in subsequent cycles.
        # For now, we just generate structures, label them, and train.
        unlabeled_ids = self.db_wrapper.get_unlabeled_ids()
        for structure_id in unlabeled_ids:
            self.label_structure_by_id(structure_id)
        self.run_training()

    def label_structure_by_id(self, structure_id: int) -> None:
        """
        Runs the labeling process for a specific structure ID.

        Args:
            structure_id: The ID of the structure to label.
        """
        print(f"Starting labeling for structure ID: {structure_id}")
        try:
            atoms = self.db_wrapper.get_atoms_by_id(structure_id)
            dft_result = self.labeling_engine.label_structure(atoms)
            self.db_wrapper.update_labels(structure_id, dft_result)
            print(f"Labeling complete for structure ID: {structure_id}")
        except Exception as e:
            import logging

            logging.error(f"Failed to label structure {structure_id}: {e}")
            self.db_wrapper.update_state(structure_id, "labeling_failed")

    def run_training(self) -> None:
        """Runs the MLIP model training process."""
        print("Starting training process.")
        try:
            training_data = self.db_wrapper.get_all_labeled_atoms()
            self.training_engine.train(training_data)
            print("Training process finished.")
        except NotImplementedError as e:
            print(f"Training not implemented: {e}")
        except Exception as e:
            import logging

            logging.error(f"Training failed: {e}")
