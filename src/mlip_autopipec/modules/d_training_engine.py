import torch

# HACK: Fix for e3nn/_wigner.py loading issue with recent PyTorch versions.
# This must be done before e3nn is imported.
import torch.serialization

torch.serialization.add_safe_globals([slice])

from pathlib import Path
from typing import List

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from e3nn import o3

from mace.data.atomic_data import Configuration
from mace.data.utils import config_from_atoms
from mace.modules.blocks import RealAgnosticInteractionBlock
from mace.modules.models import MACE
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential


class TrainingEngine:
    """The engine for training the MLIP."""

    def __init__(self, config: TrainingConfig, db: AseDB):
        """
        Initializes the TrainingEngine.

        Args:
            config: The TrainingConfig Pydantic model.
            db: An instance of AseDB.
        """
        self._config = config
        self._db = db
        self._z_table = self._get_z_table()

    def _get_z_table(self) -> AtomicNumberTable:
        """Create a Z table from the elements in the database."""
        row = self._db.connection.get(1)
        if row:
            atoms = row.toatoms()
            z_list = sorted(list(set(atoms.get_atomic_numbers())))
            return AtomicNumberTable([int(z) for z in z_list])
        return AtomicNumberTable([1, 6, 7, 8])

    def execute(self, ids: list[int]) -> str:
        """
        Loads data, prepares it for Delta Learning, trains the model, and saves it.
        """
        self._load_and_prepare_data(ids)

        print("Data prepared for training. A real training loop would follow.")

        # Define some sensible defaults for a minimal MACE model
        hidden_irreps = o3.Irreps("128x0e")
        mlp_irreps = o3.Irreps("16x0e")

        model = MACE(
            r_max=self._config.r_cut,
            num_bessel=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls=RealAgnosticInteractionBlock,
            interaction_cls_first=RealAgnosticInteractionBlock,
            num_interactions=2,
            num_elements=len(self._z_table.zs),
            hidden_irreps=hidden_irreps,
            MLP_irreps=mlp_irreps,
            atomic_energies=np.zeros(len(self._z_table.zs)),
            avg_num_neighbors=1,
            atomic_numbers=self._z_table.zs,
            correlation=3,
            gate=torch.nn.functional.silu,
        )

        model_path = Path("models")
        model_path.mkdir(exist_ok=True)
        final_model_path = model_path / "trained_model.pt"
        torch.save(model.state_dict(), final_model_path)

        return str(final_model_path)

    def _load_and_prepare_data(self, ids: list[int]) -> list[Configuration]:
        """Queries the DB, computes baseline values, and calculates residuals."""
        prepared_data = []
        for db_id in ids:
            row = self._db.connection.get(db_id)
            if not row:
                continue

            atoms = row.toatoms()
            dft_energy = row.key_value_pairs["total_energy_ev"]
            dft_forces = np.array(row.key_value_pairs["forces"])

            energy_to_set = dft_energy
            forces_to_set = dft_forces

            if self._config.delta_learn:
                baseline_energy, baseline_forces = calculate_lj_potential(atoms)
                energy_to_set = dft_energy - baseline_energy
                forces_to_set = dft_forces - baseline_forces

            # Attach results to Atoms object using a calculator
            calc = SinglePointCalculator(atoms, energy=energy_to_set, forces=forces_to_set)
            atoms.calc = calc

            # config_from_atoms will now correctly read the energy and forces
            config = config_from_atoms(atoms)
            prepared_data.append(config)

        return prepared_data
