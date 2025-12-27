# Add this to the top of the file to fix a potential pickle error with mace-torch and new torch versions
import torch.serialization
torch.serialization.add_safe_globals([slice])

import torch
import numpy as np
import os
from typing import List
from ase.calculators.singlepoint import SinglePointCalculator
from ase.atoms import Atoms

# MACE imports
from mace.modules.models import MACE
from mace.modules.blocks import RealAgnosticInteractionBlock
from mace.data.utils import config_from_atoms
from mace.tools.torch_geometric.dataloader import DataLoader
from mace.modules.loss import WeightedEnergyForcesLoss

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential


class TrainingEngine:
    def __init__(self, config: TrainingConfig, db: AseDB):
        self._config = config
        self._db = db
        self.device = "cuda" if torch.cuda.is_available() else "cpu"

    def execute(self, ids: list[int]) -> str:
        """Loads data, prepares it for Delta Learning, trains the model, and returns the path."""

        # 1. Load and prepare data
        atoms_list = self._load_and_prepare_data(ids)
        if not atoms_list:
            raise ValueError("No valid data found for the given IDs.")

        # 2. Convert to MACE format and create a DataLoader
        from mace.data.atomic_data import AtomicData
        from mace.tools import AtomicNumberTable

        atomic_numbers = np.unique(np.concatenate([atoms.get_atomic_numbers() for atoms in atoms_list])).tolist()
        z_table = AtomicNumberTable(atomic_numbers)

        configs = [config_from_atoms(atoms) for atoms in atoms_list]
        atomic_data = [AtomicData.from_config(config, z_table=z_table, cutoff=self._config.r_cut) for config in configs]
        data_loader = DataLoader(
            dataset=atomic_data,
            batch_size=len(configs), # For Cycle 1, use a single batch
        )

        # 3. Initialize MACE model with some reasonable defaults

        model_args = {
            "r_max": 5.0,
            "num_bessel": 8,
            "num_polynomial_cutoff": 5,
            "max_ell": 3,
            "interaction_cls": RealAgnosticInteractionBlock,
            "num_interactions": 2,
            "num_elements": len(atomic_numbers),
            "hidden_irreps": "128x0e + 128x1o",
            "atomic_energies": np.zeros(len(atomic_numbers)),
            "avg_num_neighbors": 8.0,
            "atomic_numbers": atomic_numbers,
            "correlation": 3,
            "gate": torch.nn.functional.silu,
        }
        model = MACE(**model_args).to(self.device)

        # 4. Set up optimizer and loss function
        optimizer = torch.optim.Adam(model.parameters(), lr=self._config.learning_rate)
        loss_fn = WeightedEnergyForcesLoss(energy_weight=1.0, forces_weight=10.0)

        # 5. Training loop
        model.train()
        for epoch in range(self._config.num_epochs):
            for batch in data_loader:
                batch = batch.to(self.device)
                optimizer.zero_grad()

                output = model(batch.to_dict())
                loss = loss_fn(output, batch)

                loss.backward()
                optimizer.step()

            print(f"Epoch {epoch+1}/{self._config.num_epochs}, Loss: {loss.item():.4f}")

        # 6. Save the model
        output_dir = "models"
        os.makedirs(output_dir, exist_ok=True)
        model_path = os.path.join(output_dir, "trained_model.pt")
        torch.save(model.state_dict(), model_path)

        return model_path

    def _load_and_prepare_data(self, ids: list[int]) -> List[Atoms]:
        """Queries DB, computes baseline values, and calculates residuals."""
        prepared_atoms_list = []
        for db_id in ids:
            atoms, kvp = self._db.get(db_id)

            if not kvp.get('was_successful', False) or atoms.calc is None:
                continue # Skip unsuccessful or unlabelled calculations

            # Retrieve data from the calculator
            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            target_energy = dft_energy
            target_forces = dft_forces

            if self._config.delta_learn:
                # For now, baseline potential is hardcoded as LJ.
                # A more advanced version would read this from config.
                base_energy, base_forces = calculate_lj_potential(atoms)
                target_energy = dft_energy - base_energy
                target_forces = dft_forces - base_forces

            # Attach results to the Atoms object via a calculator, as MACE expects
            calc = SinglePointCalculator(atoms, energy=target_energy, forces=target_forces)
            atoms.calc = calc
            prepared_atoms_list.append(atoms)

        return prepared_atoms_list
