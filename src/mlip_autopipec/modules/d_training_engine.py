# Description: The Training Engine (Module D) for MLIP model training.
import ast
from pathlib import Path
from typing import List

import numpy as np
import torch
import torch.serialization
from ase.atoms import Atoms
from mace.data import AtomicData, Configuration, config_from_atoms
from mace.modules import MACE, WeightedEnergyForcesLoss
from mace.modules.blocks import RealAgnosticInteractionBlock
from mace.tools import AtomicNumberTable, torch_geometric
from torch.optim import Adam

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import zbl_potential

# Add a safe global to handle a potential unpickling error with e3nn and torch
torch.serialization.add_safe_globals([slice])


class TrainingEngine:
    """Orchestrates the MLIP model training process."""

    def __init__(self, config: TrainingConfig, db: AseDB):
        """
        Initializes the TrainingEngine.

        Args:
            config: A TrainingConfig object with training parameters.
            db: An instance of the AseDB wrapper.
        """
        self._config = config
        self._db = db
        self.device = "cuda" if torch.cuda.is_available() else "cpu"

    def execute(self, ids: List[int]) -> str:
        """
        Trains and saves an MLIP model using the data specified by database IDs.

        Args:
            ids: A list of database IDs to use for training data.

        Returns:
            The file path to the saved, trained model.
        """
        # 1. Load and prepare data
        configs = self._load_and_prepare_data(ids)

        # If model_type is None, skip training and save a placeholder
        if self._config.model_type.lower() == "none":
            print("Model type is 'None', skipping training.")
            model_path = self._save_placeholder_model()
            return model_path

        z_table = AtomicNumberTable(
            [int(z) for z in sorted(list(set(configs[0].atomic_numbers)))]
        )
        atomic_data = [
            AtomicData.from_config(c, z_table=z_table, cutoff=self._config.r_cut)
            for c in configs
        ]

        # 2. Setup model and optimizer
        model = self._setup_model(z_table)
        loss_fn = WeightedEnergyForcesLoss(energy_weight=1.0, forces_weight=10.0)
        optimizer = Adam(model.parameters(), lr=self._config.learning_rate)

        # 3. Run training loop
        self._run_training_loop(atomic_data, model, loss_fn, optimizer)

        # 4. Save the model
        output_dir = Path("./models/")
        output_dir.mkdir(exist_ok=True)
        model_path = str(output_dir / "trained_model.pt")

        torch.save(model, model_path)

        return model_path

    def _save_placeholder_model(self) -> str:
        """Saves a simple placeholder file and returns the path."""
        output_dir = Path("./models/")
        output_dir.mkdir(exist_ok=True)
        model_path = str(output_dir / "placeholder_model.txt")
        with open(model_path, "w") as f:
            f.write("This is a placeholder model.")
        return model_path

    def _load_and_prepare_data(self, ids: List[int]) -> List[Configuration]:
        """
        Loads data from the database and applies Delta Learning adjustments.
        """
        prepared_configs = []
        for db_id in ids:
            atoms, kvp = self._db.get(db_id)
            if not kvp.get("was_successful", False):
                continue  # Skip unsuccessful calculations

            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            if self._config.delta_learn:
                baseline_energy, baseline_forces = zbl_potential(atoms)

                # The model learns the difference
                target_energy = dft_energy - baseline_energy
                target_forces = dft_forces - baseline_forces
            else:
                target_energy = dft_energy
                target_forces = dft_forces

            # Create a new Atoms object with the target values for training
            training_atoms = Atoms(
                symbols=atoms.get_chemical_symbols(),
                positions=atoms.get_positions(),
                cell=atoms.get_cell(),
                pbc=atoms.get_pbc(),
            )
            config = config_from_atoms(training_atoms)
            config.energy = target_energy
            config.forces = target_forces
            prepared_configs.append(config)

        return prepared_configs

    def _setup_model(self, z_table: AtomicNumberTable) -> MACE:
        """Initializes the MACE model with parameters from the config."""
        atomic_numbers = z_table.zs
        atomic_energies = np.zeros(len(atomic_numbers), dtype=float)
        model = MACE(
            r_max=self._config.r_cut,
            num_bessels=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls=RealAgnosticInteractionBlock,
            interaction_cls_first=RealAgnosticInteractionBlock,
            num_interactions=2,
            num_elements=len(z_table.zs),
            atomic_energies=atomic_energies,
            avg_num_neighbors=15,
            atomic_numbers=z_table.zs,
            correlation=3,
            gate=None,
            MLP_irreps="128x0e",
            hidden_irreps="128x0e",
        )
        model.to(self.device)
        return model

    def _run_training_loop(self, atomic_data, model, loss_fn, optimizer):
        """Performs the training iterations."""
        loader = torch_geometric.dataloader.DataLoader(
            dataset=atomic_data, batch_size=1, shuffle=True
        )

        for epoch in range(self._config.num_epochs):
            for batch in loader:
                optimizer.zero_grad()
                batch_data = batch.to(self.device).to_dict()

                # Forward pass
                output = model(batch_data)

                # Calculate loss
                loss = loss_fn(pred=output, ref=batch.to(self.device))

                # Backward pass and optimization
                loss.backward()
                optimizer.step()

            if epoch % 10 == 0:
                print(f"Epoch {epoch}, Loss: {loss.item()}")
