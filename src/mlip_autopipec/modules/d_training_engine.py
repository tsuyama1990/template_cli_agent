"""The Training Engine (Module D) for MLIP model training."""
import logging
from pathlib import Path
from typing import List

import torch
import numpy as np
from ase.atoms import Atoms
from mace.data import AtomicData
from mace.tools import AtomicNumberTable
from mace.modules import MACE

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils import baseline_potentials

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TrainingEngine:
    """
    Manages the training of a Machine Learning Interatomic Potential (MLIP).
    """

    def __init__(self, config: TrainingConfig, db: AseDB):
        """
        Initializes the Training Engine.

        Args:
            config: The Pydantic model containing training hyperparameters.
            db: An instance of the AseDB wrapper.
        """
        self._config = config
        self._db = db
        self._device = "cuda" if torch.cuda.is_available() else "cpu"
        logger.info(f"Using device: {self._device}")

    def execute(self, ids: List[int], model_save_path: str = "trained_model.pt") -> str:
        """
        Executes the full training workflow.

        Args:
            ids: A list of database IDs to use for training data.
            model_save_path: The path where the trained model will be saved.

        Returns:
            The path to the saved model file.
        """
        # 1. Load and prepare data
        training_data = self._load_and_prepare_data(ids)
        if not training_data:
            raise ValueError("No valid training data could be loaded.")

        # 2. Configure and train the model
        z_table = AtomicNumberTable([int(z) for z in np.unique(np.concatenate([at.get_atomic_numbers() for at in training_data]))])

        # Reverting to the full, explicit constructor as determined by introspection.
        model = MACE(
            r_max=self._config.r_cut,
            num_bessel=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls="RealAgnosticInteractionBlock", # Must be a string
            interaction_cls_first="RealAgnosticInteractionBlock", # Must be a string
            num_interactions=2,
            num_elements=len(z_table),
            hidden_irreps="128x0e + 128x1o",
            MLP_irreps="128x0e + 128x1o",
            atomic_energies=np.zeros(len(z_table)),
            avg_num_neighbors=15,
            atomic_numbers=z_table.zs,
            correlation=3,
            gate="silu",
        ).to(self._device)

        optimizer = torch.optim.Adam(model.parameters(), lr=self._config.learning_rate)

        # MACE data loader expects a list of dictionaries
        data_loader = [AtomicData.from_atoms(atoms, z_table=z_table, cutoff=self._config.r_cut).to_dict() for atoms in training_data]

        logger.info("Starting training...")
        for epoch in range(self._config.num_epochs):
            model.train()
            total_loss = 0.0
            for data in data_loader:
                optimizer.zero_grad()
                # The data dictionary needs to be moved to the correct device
                data_on_device = {k: v.to(self._device) if hasattr(v, 'to') else v for k, v in data.items()}
                output = model(data_on_device)

                # Ensure target tensors are also on the correct device
                energy_target = data_on_device["energy"]
                forces_target = data_on_device["forces"]

                loss = (output["energy"] - energy_target).pow(2).mean() + (output["forces"] - forces_target).pow(2).mean()
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

            avg_loss = total_loss / len(data_loader)
            logger.info(f"Epoch {epoch+1}/{self._config.num_epochs}, Average Loss: {avg_loss:.6f}")

        # 3. Save the model
        Path(model_save_path).parent.mkdir(parents=True, exist_ok=True)
        torch.save(model.state_dict(), model_save_path)
        logger.info(f"Model saved to {model_save_path}")

        return model_save_path

    def _load_and_prepare_data(self, ids: List[int]) -> List[Atoms]:
        """
        Loads data from the DB and prepares it for Delta Learning.
        """
        prepared_atoms_list = []
        for db_id in ids:
            atoms, dft_result = self._db.get(db_id)
            if not dft_result.was_successful:
                logger.warning(f"Skipping unsuccessful calculation with ID {db_id}.")
                continue

            target_energy = dft_result.total_energy_ev
            target_forces = np.array(dft_result.forces)

            if self._config.delta_learn:
                baseline_energy, baseline_forces = baseline_potentials.calculate_lennard_jones(atoms)
                target_energy -= baseline_energy
                target_forces -= baseline_forces

            # MACE expects the target values to be in the .info dictionary as tensors
            atoms.info["energy"] = torch.tensor(target_energy, dtype=torch.float64)
            atoms.info["forces"] = torch.tensor(target_forces, dtype=torch.float64)
            prepared_atoms_list.append(atoms)

        return prepared_atoms_list
