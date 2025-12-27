"""The Training Engine for MLIP models."""
import os
import torch
from typing import List
from ..data.database import AseDB
from ..data.models import TrainingConfig

class TrainingEngine:
    """A dummy Training Engine that creates a placeholder model file."""

    def __init__(self, config: TrainingConfig, db: AseDB, model_dir: str = "models"):
        self._config = config
        self._db = db
        self._model_dir = model_dir
        os.makedirs(self._model_dir, exist_ok=True)

    def execute(self, ids: List[int]) -> str:
        """
        Creates a dummy model file.
        """
        # Ensure at least one valid ID was passed
        if not any(self._db.get(db_id) for db_id in ids):
            raise ValueError("No valid database IDs provided for training.")

        model_path = os.path.join(self._model_dir, "dummy_model.pt")

        # Create a simple dummy tensor and save it
        dummy_tensor = torch.tensor([1, 2, 3])
        torch.save(dummy_tensor, model_path)

        return model_path
