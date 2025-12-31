import os
import torch
from ase import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.labelling_engine import LabellingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine


class Orchestrator:
    """
    Contains the main business logic for the MLIP-AutoPipe workflow.
    """

    def __init__(self, db_path: str, qe_command: str, pseudo_dir: str):
        self.db = AseDB(db_path)
        self.labelling_engine = LabellingEngine(qe_command, pseudo_dir)
        self.training_engine = TrainingEngine()

    def run_cycle01_workflow(self, atoms: Atoms):
        """
        Executes the simple, linear workflow for Cycle 1.
        """
        print("Starting Cycle 01 Workflow...")

        # 1. Add structure to database
        atoms_id = self.db.add_atoms(atoms, state="awaiting_labelling")
        print(f"Added structure {atoms_id} to the database.")

        # 2. Run Labelling Engine
        print("Running Labelling Engine...")
        dft_result = self.labelling_engine.run(atoms)

        if dft_result.status == "failed":
            self.db.update_state(atoms_id, "labelling_failed")
            print(f"Labelling failed for structure {atoms_id}.")
            return

        # 3. Write results to database
        self.db.write_dft_result(atoms_id, dft_result)
        print("DFT calculation successful. Results written to database.")

        # 4. Get training data
        training_data = self.db.get_training_data()
        if not training_data[0]:
            print("No training data found. Aborting.")
            return

        # 5. Run Training Engine
        print("Running Training Engine...")
        model = self.training_engine.train(training_data)
        print("Model training complete.")

        # 6. Save model
        model_path = "model_cycle01.pt"
        torch.save(model.state_dict(), model_path)
        print(f"Model saved to {model_path}")
