import os

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig


class TrainingEngine:
    """
    (Placeholder) Manages the MLIP training workflow.

    This is a placeholder implementation that does not perform any actual training.
    It simulates a successful run by creating an empty model file. The real
    MACE integration has been deferred to a future cycle due to API instability.
    """
    def __init__(self, config: TrainingConfig, db: AseDB):
        self._config = config
        self._db = db
        print("Placeholder TrainingEngine initialized.")

    def execute(self, ids: list[int], output_dir: str = "models") -> str:
        """
        Simulates a training run by creating an empty model file.
        """
        print("Executing placeholder training...")

        has_valid_data = any(
            self._db.get(db_id)[1].get("was_successful", False) for db_id in ids
        )

        if not has_valid_data:
            raise ValueError(
                "No successful DFT calculations found in the given IDs. "
                "Skipping training."
            )

        os.makedirs(output_dir, exist_ok=True)
        model_path = os.path.join(output_dir, "mace_model_placeholder.pt")

        with open(model_path, "w") as f:
            f.write("This is a placeholder model.\n")

        print(f"Placeholder training complete. Model saved to: {model_path}")
        return model_path
