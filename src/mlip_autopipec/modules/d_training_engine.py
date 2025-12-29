
import torch
from mace.data.atomic_data import AtomicData
from mace.data.utils import config_from_atoms
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils import baseline_potentials


class TrainingEngine:
    def __init__(self, config: TrainingConfig, db: AseDB):
        self._config = config
        self._db = db
        self._z_table = None

    def execute(self, ids: list[int]) -> str:
        """
        Loads data for given IDs from the DB, prepares it for Delta Learning,
        trains the model, and returns the path to the saved model file.
        """
        prepared_data = self._load_and_prepare_data(ids)
        # In a real scenario, a full training loop would be here.
        # We mock the actual training for this cycle's tests.
        model_path = self._save_model(prepared_data)
        return model_path

    def _load_and_prepare_data(self, ids: list[int]) -> list[AtomicData]:
        """
        Queries the DB, computes baseline values if delta learning is enabled,
        and prepares the data in the format expected by MACE.
        """
        all_atoms = []
        for db_id in ids:
            atoms, dft_kvp = self._db.get(db_id)
            if not dft_kvp.get("was_successful", False):
                continue

            # The energy/forces are stored on the Atoms object by the calculator
            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            if self._config.delta_learn:
                # For now, hardcode LJ params. In the future, this would come
                # from the Heuristic Engine.
                epsilon = 0.0103  # eV for Ar
                sigma = 3.4  # Angstrom for Ar

                baseline_energy, baseline_forces = baseline_potentials.calculate_lj(
                    atoms, epsilon=epsilon, sigma=sigma
                )
                # Store the final delta values in info/arrays
                atoms.info["energy"] = dft_energy - baseline_energy
                atoms.arrays["forces"] = dft_forces - baseline_forces
            else:
                atoms.info["energy"] = dft_energy
                atoms.arrays["forces"] = dft_forces

            all_atoms.append(atoms)

        # Convert ASE Atoms to MACE AtomicData
        if not self._z_table:
            all_zs = set()
            for atoms in all_atoms:
                all_zs.update(atoms.get_atomic_numbers())
            self._z_table = AtomicNumberTable([int(z) for z in sorted(list(all_zs))])

        atomic_data_list = []
        for atoms in all_atoms:
            config = config_from_atoms(atoms)
            atomic_data = AtomicData.from_config(
                config, z_table=self._z_table, cutoff=self._config.r_cut
            )
            # Manually set the final energy and forces for MACE
            atomic_data.energy = torch.tensor(atoms.info["energy"], dtype=torch.float64)
            atomic_data.forces = torch.tensor(atoms.arrays["forces"], dtype=torch.float64)
            atomic_data_list.append(atomic_data)

        return atomic_data_list

    def _save_model(self, data) -> str:
        # Placeholder for the actual training and saving logic
        model_path = "models/trained_model.pt"
        # Ensure the directory exists
        import os
        os.makedirs(os.path.dirname(model_path), exist_ok=True)
        print(f"Pretending to train and save model to {model_path}")
        return model_path
