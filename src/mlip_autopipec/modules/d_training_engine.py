import torch
import os
import numpy as np
from typing import List
from ase.db import connect

from mace.tools import AtomicNumberTable
from mace.modules.models import MACE

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential, calculate_zbl_potential
from mlip_autopipec.utils.mace_utils import to_atomic_data

class TrainingEngine:
    def __init__(self, config: TrainingConfig, db: AseDB):
        self.config = config
        self.db = db
        self.z_table = None

    def _load_and_prepare_data(self, ids: List[int]) -> List: # Returns List[AtomicData]
        """
        Loads data from the database, calculates baseline potentials,
        and prepares it in the format required by MACE.
        """
        all_atoms = []
        with connect(self.db.db_path) as db:
            for db_id in ids:
                row = db.get(id=db_id)
                atoms = row.toatoms()
                dft_energy = row.energy
                dft_forces = row.forces

                if self.config.delta_learn:
                    if self.config.baseline_potential == 'lj':
                        baseline_energy, baseline_forces = calculate_lj_potential(atoms)
                    elif self.config.baseline_potential == 'zbl':
                        baseline_energy, baseline_forces = calculate_zbl_potential(atoms)
                    else:
                        baseline_energy, baseline_forces = 0.0, np.zeros((len(atoms), 3))

                    atoms.calc = None # Use modern syntax
                    atoms.info['energy'] = dft_energy - baseline_energy
                    atoms.arrays['forces'] = dft_forces - baseline_forces
                else:
                    atoms.calc = None # Use modern syntax
                    atoms.info['energy'] = dft_energy
                    atoms.arrays['forces'] = dft_forces

                all_atoms.append(atoms)

        atomic_numbers = set()
        for atoms in all_atoms:
            atomic_numbers.update(atoms.get_atomic_numbers())

        self.z_table = AtomicNumberTable([int(z) for z in sorted(list(atomic_numbers))])

        atomic_data_list = [
            to_atomic_data(atoms, z_table=self.z_table, cutoff=self.config.r_cut)
            for atoms in all_atoms
        ]

        return atomic_data_list

    def execute(self, ids: List[int], output_dir: str = "models") -> str:
        """
        Loads data, trains the MACE model, and saves it.
        Returns the path to the saved model file.
        """

        training_data = self._load_and_prepare_data(ids)

        model = MACE(
            r_max=self.config.r_cut,
            num_bessels=8,
            num_polynomial_basis=4,
            atomic_numbers=self.z_table.zs,
            correlation=3,
            hidden_irreps='16x0e + 16x1o',
            MLP_irreps='16x0e',
            atomic_energies=np.zeros(len(self.z_table.zs)),
            avg_num_neighbors=15.0,
            atomic_inter_scale=1.0,
            atomic_inter_shift=0.0
        )

        optimizer = torch.optim.Adam(model.parameters(), lr=self.config.learning_rate)

        for epoch in range(self.config.num_epochs):
            for data in training_data:
                optimizer.zero_grad()
                output = model(data.to_dict())
                loss = torch.sum((output['energy'] - data.energy)**2)
                loss.backward()
                optimizer.step()
            print(f"Epoch {epoch+1}/{self.config.num_epochs}, Loss: {loss.item()}")

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        model_path = os.path.join(output_dir, "trained_model.pt")
        torch.save(model, model_path)

        return model_path
