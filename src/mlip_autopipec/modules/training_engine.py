"""The Training Engine for MLIP models."""

import numpy as np
import torch
from ase import Atoms
from ase.calculators.lj import LennardJones
from mace.data.atomic_data import AtomicData
from mace.data.utils import config_from_atoms
from mace.modules.blocks import RealAgnosticInteractionBlock
from mace.modules.models import MACE
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.database import AseDB


class TrainingEngine:
    """Manages the training of the MLIP model."""

    def __init__(self, db: AseDB, config: dict):
        """
        Initializes the TrainingEngine.
        Args:
            db: An AseDB instance for database interaction.
            config: A dictionary containing configuration for the engine.
        """
        self.db = db
        self.config = config
        self.z_table = None

    def run(self):
        """
        Runs the training process.
        """
        print("Starting training process...")
        configs, dft_results = self._prepare_dataset()
        if not configs:
            print("No labelled data found. Skipping training.")
            return

        baseline_results = self._calculate_baseline(configs)
        delta_results = self._calculate_delta(dft_results, baseline_results)
        self._train_model(configs, delta_results)
        print("Training process completed.")

    def _prepare_dataset(self) -> (list, list):
        """
        Fetches labelled data from the database and prepares it for training.
        Returns a tuple of (configurations, dft_results).
        """
        labelled_atoms = self.db.get_atoms_by_state("labelled")

        atomic_numbers = set()
        for atoms in labelled_atoms:
            atomic_numbers.update(atoms.get_atomic_numbers())

        if not atomic_numbers:
            return [], []

        self.z_table = AtomicNumberTable([int(z) for z in sorted(list(atomic_numbers))])

        configs = [config_from_atoms(atoms) for atoms in labelled_atoms]
        dft_results = [
            {"energy": atoms.get_potential_energy(), "forces": atoms.get_forces()}
            for atoms in labelled_atoms
        ]
        return configs, dft_results

    def _calculate_baseline(self, configs: list) -> list:
        """Calculates the baseline energy and forces using Lennard-Jones."""
        lj = LennardJones(sigma=3.4, epsilon=0.0103)
        baseline_results = []
        for config in configs:
            atoms = Atoms(
                numbers=config.atomic_numbers,
                positions=config.positions,
                cell=config.cell,
                pbc=config.pbc,
            )
            atoms.calc = lj
            baseline_results.append(
                {"energy": atoms.get_potential_energy(), "forces": atoms.get_forces()}
            )
        return baseline_results

    def _calculate_delta(self, dft_results: list, baseline_results: list) -> list:
        """Calculates the delta between DFT and baseline."""
        delta_results = []
        for dft, baseline in zip(dft_results, baseline_results, strict=True):
            delta_results.append(
                {
                    "energy": dft["energy"] - baseline["energy"],
                    "forces": dft["forces"] - baseline["forces"],
                }
            )
        return delta_results

    def _train_model(self, configs: list, delta_results: list):
        """
        Trains the MACE model on the delta dataset.
        """
        cutoff = self.config.get("cutoff", 5.0)
        atomic_data = [
            AtomicData.from_config(config, z_table=self.z_table, cutoff=cutoff)
            for config in configs
        ]

        for ad, delta in zip(atomic_data, delta_results, strict=True):
            ad.energy = torch.tensor(delta["energy"], dtype=torch.float32)
            ad.forces = torch.tensor(delta["forces"], dtype=torch.float32)

        model_config = {
            "r_max": cutoff,
            "num_bessel": 8,
            "num_polynomial_cutoff": 5,
            "max_ell": 3,
            "interaction_cls": RealAgnosticInteractionBlock,
            "num_interactions": 2,
            "hidden_irreps": "128x0e + 128x1o",
            "atomic_energies": np.zeros(len(self.z_table)),
            "avg_num_neighbors": 15.0,
            "atomic_numbers": self.z_table.zs,
            "correlation": 3,
        }

        model = MACE(**model_config)

        print("Skipping actual training loop for Cycle 1. Saving untrained model.")

        output_path = self.config.get("model_path", "model.pt")
        torch.save(model.state_dict(), output_path)
        print(f"Model saved to {output_path}")
