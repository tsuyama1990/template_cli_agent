from typing import List
import json
import numpy as np
from ase.db.row import AtomsRow
from pathlib import Path

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig, DFTResult
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential

class TrainingEngine:
    """
    Manages the MLIP model training workflow.
    """
    def __init__(self, config: TrainingConfig, db: AseDB):
        """
        Initializes the TrainingEngine.

        Args:
            config: The Pydantic model containing training hyperparameters.
            db: An instance of the AseDB class for data access.
        """
        self._config = config
        self._db = db

    def execute(self, ids: List[int]) -> str:
        """
        Loads data, prepares it for Delta Learning, trains the model,
        and returns the path to the saved model file.

        Args:
            ids: A list of database IDs to use for training.

        Returns:
            The path to the saved model file.
        """
        prepared_data = self._load_and_prepare_data(ids)

        print(f"Starting training with {len(prepared_data)} data points...")
        print(f"Config: {self._config}")

        # Mock training loop
        for epoch in range(self._config.num_epochs):
            pass

        # Mock model saving: create the directory and a dummy file.
        model_path_str = f"models/{self._config.model_type}_model.pt"
        model_path = Path(model_path_str)
        model_path.parent.mkdir(parents=True, exist_ok=True)
        with open(model_path, "w") as f:
            f.write("This is a mock model file.")

        print(f"Training complete. Model saved to: {model_path_str}")

        return model_path_str

    def _load_and_prepare_data(self, ids: List[int]) -> List[dict]:
        """
        Queries the database, computes baseline values if Delta Learning is enabled,
        and calculates the residuals for the ML model to learn.
        """
        prepared_data = []
        for db_id in ids:
            row = self._db.get(db_id)
            if not row or not row.key_value_pairs["was_successful"]:
                continue

            atoms = row.toatoms()
            dft_result = DFTResult.model_validate(row.data)

            target_energy = dft_result.total_energy_ev
            target_forces = np.array(dft_result.forces)

            if self._config.delta_learn:
                if self._config.baseline_potential == 'lj':
                    baseline_energy, baseline_forces = calculate_lj_potential(atoms)
                else:
                    raise ValueError(f"Unknown baseline potential: {self._config.baseline_potential}")

                target_energy -= baseline_energy
                target_forces -= baseline_forces

            prepared_data.append({
                "atoms": atoms,
                "target_energy": target_energy,
                "target_forces": target_forces,
            })

        return prepared_data
