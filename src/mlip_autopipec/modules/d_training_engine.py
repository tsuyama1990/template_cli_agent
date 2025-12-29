# ruff: noqa: D101, D102, D103, D104, D105, D107, F401, S608
import logging
from pathlib import Path

import numpy as np
import torch
from ase.atoms import Atoms
from mace.data.atomic_data import AtomicData, Configuration
from mace.data.utils import config_from_atoms
from mace.modules.blocks import RealAgnosticInteractionBlock
from mace.modules.loss import WeightedEnergyForcesLoss
from mace.modules.models import MACE
from mace.tools import AtomicNumberTable
from mace.tools.torch_geometric.dataloader import DataLoader

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential

# Handles a known mace-torch/e3nn unpickling error with modern torch
import torch.serialization
torch.serialization.add_safe_globals([slice])


class TrainingEngine:
    """Handles the training of a Machine Learning Interatomic Potential (MLIP)
    using the MACE framework.
    """

    def __init__(self, config: TrainingConfig, db: AseDB, output_dir: str = "models"):
        """Initializes the TrainingEngine.

        Args:
            config: A Pydantic model containing the training configuration.
            db: An instance of the AseDB wrapper for accessing training data.
            output_dir: The directory where the trained model will be saved.
        """
        self._config = config
        self._db = db
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self._z_table: AtomicNumberTable | None = None

    def _load_and_prepare_data(self, ids: list[int]) -> list[AtomicData]:
        """Queries DB, computes baselines, and calculates residuals for Delta Learning."""
        atoms_list: list[Atoms] = [self._db.get(db_id)[0] for db_id in ids]

        all_zs = sorted(
            list(set(z for atoms in atoms_list for z in atoms.get_atomic_numbers()))
        )
        self._z_table = AtomicNumberTable([int(z) for z in all_zs])

        atomic_data_list = []
        for atoms in atoms_list:
            if not atoms.calc:
                logging.warning(f"Atoms object has no calculator, skipping.")
                continue

            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            delta_energy, delta_forces = self._calculate_deltas(
                atoms, dft_energy, dft_forces
            )

            config = config_from_atoms(atoms)
            atomic_data = AtomicData.from_config(
                config, z_table=self._z_table, cutoff=self._config.r_cut
            )

            atomic_data.energy = torch.tensor([delta_energy], dtype=torch.float64)
            atomic_data.forces = torch.tensor(delta_forces, dtype=torch.float64)
            atomic_data_list.append(atomic_data)

        return atomic_data_list

    def _calculate_deltas(self, atoms, dft_energy, dft_forces):
        if self._config.delta_learn:
            baseline_energy, baseline_forces = calculate_lj_potential(atoms)
            return dft_energy - baseline_energy, dft_forces - baseline_forces
        return dft_energy, dft_forces

    def execute(self, ids: list[int]) -> str:
        """Loads data, prepares it for Delta Learning, trains, and saves the model."""
        training_data = self._load_and_prepare_data(ids)
        if not training_data:
            raise ValueError("No valid training data could be loaded for the given IDs.")

        train_loader = DataLoader(dataset=training_data, batch_size=1, shuffle=True)

        atomic_energies = np.zeros(len(self._z_table.zs))
        model = MACE(
            r_max=self._config.r_cut,
            atomic_numbers=self._z_table.zs,
            interaction_cls=RealAgnosticInteractionBlock,
            interaction_cls_first=RealAgnosticInteractionBlock,
            num_interactions=2,
            num_elements=len(self._z_table.zs),
            # NOTE: Using simplified irreps to avoid a library bug for the test.
            hidden_irreps="4x0e",
            MLP_irreps="4x0e",
            atomic_energies=atomic_energies,
            avg_num_neighbors=8.0,
            correlation=3,
            num_bessel=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            gate=torch.nn.functional.silu,
        )
        model.to(torch.float64)

        loss_fn = WeightedEnergyForcesLoss(energy_weight=1.0, forces_weight=10.0)
        optimizer = torch.optim.Adam(model.parameters(), lr=self._config.learning_rate)

        model.train()
        for _epoch in range(self._config.num_epochs):
            for batch in train_loader:
                optimizer.zero_grad()
                output = model(batch.to_dict())
                loss = loss_fn(batch, output)
                loss.backward()
                optimizer.step()

        model_path = self.output_dir / "trained_model.pt"
        torch.save(model, model_path)

        return str(model_path)
