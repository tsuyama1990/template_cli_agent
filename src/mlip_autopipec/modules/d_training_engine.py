import subprocess
import shutil
from pathlib import Path
from typing import List, Optional
import logging

import numpy as np
from ase.io import write

from ..data.database import AseDB
from ..data.models import TrainingConfig
from ..utils.baseline_potentials import calculate_lj_potential

class TrainingEngine:
    def __init__(self, config: TrainingConfig, db: AseDB):
        self._config = config
        self._db = db

    def execute(self, ids: List[int], model_save_path: str = "trained_model.pt", work_dir: Path = Path("./mace_temp")) -> Optional[str]:
        if not ids:
            logging.error("Cannot train on an empty list of database IDs.")
            return None

        work_dir.mkdir(exist_ok=True)
        train_file = work_dir / "train.xyz"

        try:
            training_atoms = self._prepare_data_for_xyz(ids)

            if not training_atoms:
                logging.error("No valid training data found for the given IDs.")
                return None

            write(train_file, training_atoms)

            command = self._build_cli_command(train_file, model_save_path)

            logging.info(f"Running MACE training command: {' '.join(command)}")

            process_result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                cwd=work_dir
            )

            if process_result.returncode != 0:
                logging.error(f"MACE training failed. Error:\n{process_result.stderr}")
                return None

            logging.info(f"Training complete. Model saved to: {model_save_path}")
            return model_save_path

        finally:
            if work_dir.exists():
                shutil.rmtree(work_dir)

    def _prepare_data_for_xyz(self, ids: List[int]) -> List:
        prepared_atoms_list = []
        for db_id in ids:
            entry = self._db.get(db_id)
            if not entry or not entry[1].was_successful:
                continue

            atoms, dft_result = entry
            if dft_result.total_energy_ev is None or dft_result.forces is None:
                continue

            atoms_copy = atoms.copy()
            if self._config.delta_learn:
                baseline_energy, baseline_forces = calculate_lj_potential(atoms_copy)
                energy_target = dft_result.total_energy_ev - baseline_energy
                forces_target = np.array(dft_result.forces) - baseline_forces
            else:
                energy_target = dft_result.total_energy_ev
                forces_target = np.array(dft_result.forces)

            atoms_copy.info['energy'] = energy_target
            atoms_copy.arrays['forces'] = forces_target
            atoms_copy.calc = None

            prepared_atoms_list.append(atoms_copy)

        return prepared_atoms_list

    def _build_cli_command(self, train_file: Path, model_save_path: str) -> List[str]:
        command = [
            "mace_run_train",
            f"--name=MLIP_AutoPipe_MACE",
            f"--train_file={train_file}",
            f"--model_dir={Path(model_save_path).parent}",
            f"--hidden_irreps=128x0e+128x1o",
            f"--r_max={self._config.r_cut}",
            f"--num_epochs={self._config.num_epochs}",
            f"--learning_rate={self._config.learning_rate}",
        ]
        return command
