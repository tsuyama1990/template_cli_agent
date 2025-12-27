# src/mlip_autopipec/modules/d_training_engine.py
"""
This module defines the TrainingEngine, which is responsible for training
the Machine Learning Interatomic Potential (MLIP).
"""
import subprocess
from pathlib import Path
from typing import List

from ase.atoms import Atoms
from ase.io import write
import numpy as np

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils import baseline_potentials


class TrainingEngine:
    """
    This class encapsulates the logic for training a Machine Learning Interatomic
    Potential (MLIP) model. It handles data loading, preparation for Delta
    Learning, and interfacing with the underlying MLIP framework (e.g., MACE)
    via its command-line interface.
    """

    def __init__(self, config: TrainingConfig, db: AseDB):
        self._config = config
        self._db = db
        self._work_dir = Path.cwd()
        self._model_name = "trained_model"
        self._model_path = self._work_dir / f"{self._model_name}.model"

    def execute(self, ids: List[int]) -> Path:
        """
        Executes the full training workflow.

        Args:
            ids: A list of database IDs to be used as the training set.

        Returns:
            The file path to the saved, trained model.
        """
        training_data = self._load_and_prepare_data(ids)
        if not training_data:
            raise ValueError("No valid training data found for the given IDs.")

        train_file_path = self._work_dir / "temp_train.xyz"
        try:
            write(train_file_path, training_data)

            command = [
                "mace_run_train",
                f"--name={self._model_name}",
                f"--train_file={train_file_path}",
                "--energy_key=energy",
                "--forces_key=forces",
                f"--r_max={self._config.r_cut}",
                f"--max_num_epochs={self._config.num_epochs}",
                f"--lr={self._config.learning_rate}",
                "--hidden_irreps=128x0e",
                "--device=cpu",
                "--save_cpu",
                "--E0s=average",
                "--batch_size=1",
                "--valid_batch_size=1",
            ]

            result = subprocess.run(command, capture_output=True, text=True, check=False)

            if result.returncode != 0:
                error_message = (
                    f"MACE training failed with exit code {result.returncode}.\\n"
                    f"Stdout:\\n{result.stdout}\\n"
                    f"Stderr:\\n{result.stderr}"
                )
                raise RuntimeError(error_message)

        finally:
            if train_file_path.exists():
                train_file_path.unlink()

        final_model_path = self._work_dir / f"{self._model_name}.model"
        final_pt_path = self._work_dir / "models" / "trained_model.pt"
        final_pt_path.parent.mkdir(exist_ok=True)
        if final_model_path.exists():
            final_model_path.rename(final_pt_path)
            return final_pt_path

        raise FileNotFoundError("MACE did not produce the expected model file.")

    def _load_and_prepare_data(self, ids: List[int]) -> List[Atoms]:
        """
        Loads data from the database and prepares it for training.
        """
        prepared_atoms_list = []
        for db_id in ids:
            row = self._db.get(db_id)
            if not row or not row.get("was_successful"):
                continue

            atoms = Atoms(
                numbers=row["numbers"],
                positions=row["positions"],
                cell=row["cell"],
                pbc=row["pbc"],
            )
            dft_energy = row["total_energy_ev"]
            dft_forces = np.array(row["forces"])

            target_energy = dft_energy
            target_forces = dft_forces

            if self._config.delta_learn:
                if self._config.baseline_potential == "lj":
                    baseline_energy, baseline_forces = baseline_potentials.calculate_lj_potential(atoms)
                else:
                    raise ValueError("Unsupported baseline potential specified.")

                target_energy -= baseline_energy
                target_forces -= baseline_forces

            atoms.info["energy"] = target_energy
            atoms.arrays["forces"] = target_forces
            prepared_atoms_list.append(atoms)

        return prepared_atoms_list
