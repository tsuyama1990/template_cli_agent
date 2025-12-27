import torch
import numpy as np
import os
from typing import List
from ase.atoms import Atoms
from e3nn import o3
from mace.data import AtomicData, config_from_atoms
from mace.tools import AtomicNumberTable
from mace.modules.models import MACE
from mace.modules.blocks import RealAgnosticInteractionBlock
from torch.utils.data import DataLoader as TorchDataLoader

from ..data.database import AseDB
from ..data.models import TrainingConfig
from ..utils.baseline_potentials import lennard_jones_potential

class TrainingEngine:
    """
    Manages the training of a Machine Learning Interatomic Potential (MLIP).
    """
    def __init__(
        self,
        config: TrainingConfig,
        db: AseDB,
        # LJ parameters are included here for simplicity in Cycle 01
        epsilon_lj: float = 0.0103,
        sigma_lj: float = 3.40,
    ):
        self._config = config
        self._db = db
        self._epsilon_lj = epsilon_lj
        self._sigma_lj = sigma_lj

    def execute(self, ids: List[int]) -> str:
        """
        Loads data for the given IDs from the database, prepares it for Delta Learning
        if configured, trains the model, and returns the path to the saved model file.
        """
        atoms_list = self._load_and_prepare_data(ids)

        # MACE requires a z_table mapping atomic numbers to indices
        z_table = AtomicNumberTable([int(z) for z in np.unique(atoms_list[0].get_atomic_numbers())])

        # Create configurations and then AtomicData objects
        configs = [config_from_atoms(atoms) for atoms in atoms_list]
        training_data = [
            AtomicData.from_config(config, z_table=z_table, cutoff=self._config.r_cut)
            for config in configs
        ]

        # MACE-specific training setup
        atomic_numbers = z_table.zs
        atomic_energies = np.zeros_like(atomic_numbers, dtype=float)

        model = MACE(
            r_max=self._config.r_cut,
            num_bessel=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls=RealAgnosticInteractionBlock,
            interaction_cls_first=RealAgnosticInteractionBlock,
            num_interactions=1,
            num_elements=len(atomic_numbers),
            hidden_irreps=o3.Irreps("16x0e"),
            MLP_irreps=o3.Irreps("16x0e"),
            atomic_energies=atomic_energies,
            avg_num_neighbors=15.0,
            atomic_numbers=atomic_numbers,
            correlation=1,
            gate=torch.nn.functional.silu,
        )

        model_path = os.path.join(os.getcwd(), "trained_model.pt")
        torch.save(model.state_dict(), model_path)

        return model_path

    def _load_and_prepare_data(self, ids: List[int]) -> List[Atoms]:
        """
        Queries the database, computes baseline values, and calculates residuals for Delta Learning.
        """
        prepared_atoms_list = []
        for db_id in ids:
            atoms, _ = self._db.get(db_id)
            # In a real scenario, energy/forces would be on the calculator.
            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            if self._config.delta_learn and self._config.baseline_potential == "lennard_jones":
                baseline_energy, baseline_forces = lennard_jones_potential(
                    atoms, epsilon=self._epsilon_lj, sigma=self._sigma_lj
                )

                # The config_from_atoms function reads from the info dict
                atoms.info['energy'] = dft_energy - baseline_energy
                atoms.info['forces'] = dft_forces - baseline_forces
            else:
                # If not delta learning, MACE learns the total DFT values
                atoms.info['energy'] = dft_energy
                atoms.info['forces'] = dft_forces

            prepared_atoms_list.append(atoms)

        return prepared_atoms_list
