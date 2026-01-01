# src/mlip_autopipec/modules/d_training_engine.py

from pathlib import Path

import ase.calculators.lj
import numpy as np
from ase import Atoms

from mlip_autopipec.configs.models import MLIPTrainingConfig
from mlip_autopipec.data.database import AseDBWrapper


class TrainingEngine:
    """An engine for training an MLIP model."""

    def __init__(self, config: MLIPTrainingConfig, db_wrapper: AseDBWrapper):
        """
        Initializes the TrainingEngine.

        Args:
            config: The MLIP training configuration.
            db_wrapper: The database wrapper instance.
        """
        self.config = config
        self.db_wrapper = db_wrapper

    def _calculate_baseline(self, atoms: Atoms) -> dict:
        """Calculates energy and forces using a baseline potential."""
        calculator = ase.calculators.lj.LennardJones()
        atoms.calc = calculator
        return {
            "energy": atoms.get_potential_energy(),
            "forces": atoms.get_forces(),
        }

    def _run_actual_training(self, training_data: list[Atoms]):
        """(Mocked) Runs the actual training library and saves the model."""
        # --- Mocked Training Call ---
        # In a real implementation, this is where you would interface with a library
        # like MACE or ACE.
        # e.g., ace_trainer.train(training_data)
        print("Mocking call to MLIP training library...")
        print("Training complete.")
        # --- End Mocked Training Call ---

        # Save the final model artifact
        model_path = Path(f"model.{self.config.model_type}")
        model_path.touch()
        print(f"Saved trained model to {model_path}")

    def run(self):
        """Runs the training process."""
        labeled_atoms_data = self.db_wrapper.get_labeled_atoms()
        print(f"Found {len(labeled_atoms_data)} labeled structures for training.")

        training_data = []
        for atoms, kvp in labeled_atoms_data:
            # The data from the db is often a list, convert to numpy array
            dft_forces = np.array(kvp.get("data", {}).get("forces"))
            dft_energy = kvp.get("data", {}).get("energy")

            if dft_energy is None or dft_forces.size == 0:
                print(f"Skipping structure {kvp.get('id')} due to missing data.")
                continue

            processed_atoms = atoms.copy()
            processed_atoms.info["dft_energy"] = dft_energy
            processed_atoms.info["dft_forces"] = dft_forces

            if self.config.delta_learning:
                baseline_results = self._calculate_baseline(processed_atoms.copy())
                target_energy = dft_energy - baseline_results["energy"]
                target_forces = dft_forces - baseline_results["forces"]
            else:
                target_energy = dft_energy
                target_forces = dft_forces

            processed_atoms.info["target_energy"] = target_energy
            processed_atoms.info["target_forces"] = target_forces
            training_data.append(processed_atoms)

        if not training_data:
            print("No valid data available for training. Exiting.")
            return

        print(f"Prepared {len(training_data)} structures for the training library.")
        self._run_actual_training(training_data)
