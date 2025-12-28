import os
import subprocess
import tempfile

from ase.io import write as ase_write

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig


class TrainingEngine:
    """
    The Training Engine (Module D) for MLIP model training, using the `mace_run_train` CLI.
    """
    def __init__(self, config: TrainingConfig, db: AseDB):
        """
        Initializes the TrainingEngine.
        """
        self._config = config
        self._db = db

    def execute(self, ids: list[int], model_save_path: str = "models/trained_model.pt") -> str:
        """
        Loads data, trains the MLIP model using mace_run_train, and saves it.
        """
        if not ids:
            raise ValueError("No database IDs provided for training.")

        with tempfile.TemporaryDirectory() as tmpdir:
            train_file_path = os.path.join(tmpdir, "train.xyz")

            # 1. Prepare training data
            atoms_list = [self._db.get(i)[0] for i in ids if self._db.get(i)[0] is not None]
            if not atoms_list:
                raise ValueError("No valid atoms objects found for the given IDs.")

            ase_write(train_file_path, atoms_list, format="extxyz", append=False)

            # 2. Construct mace_run_train command
            command = [
                "mace_run_train",
                f"--name={os.path.basename(model_save_path)}",
                f"--train_file={train_file_path}",
                "--model=MACE",
                "--hidden_irreps='128x0e + 128x1o'",
                f"--r_max={self._config.r_cut}",
                "--batch_size=5", # Keep it small for demonstration
                f"--max_num_epochs={self._config.num_epochs}",
                "--device=cpu", # Assuming CPU for simplicity
                "--default_dtype=float64",
                f"--results_dir={os.path.dirname(model_save_path)}",
            ]

            # 3. Execute training
            try:
                subprocess.run(
                    command,
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"MACE training failed with exit code {e.returncode}.\n"
                    f"Stderr: {e.stderr}\nStdout: {e.stdout}"
                ) from e
            except FileNotFoundError as e:
                 raise RuntimeError(
                    "mace_run_train command not found. "
                    "Is the mace-torch package installed and in your PATH?"
                 ) from e

        return model_save_path
