import numpy as np
import torch
from ase.atoms import Atoms
from mace.data.atomic_data import AtomicData
from mace.data.utils import config_from_atoms
from mace.modules.models import MACE
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import MLIPTraining


class TrainingEngine:
    """Engine for training an MLIP model."""

    def __init__(self, db_wrapper: AseDBWrapper, config: MLIPTraining):
        self.db_wrapper = db_wrapper
        self.config = config
        self.z_table = None

    def execute(self):
        """Executes the training workflow."""
        labeled_rows = self.db_wrapper.connection.select(labeled=True)

        atoms_list = [row.toatoms() for row in labeled_rows]

        if not atoms_list:
            print("No labeled structures found to train on.")
            return

        self.z_table = self._create_z_table(atoms_list)

        train_data = self._prepare_data(atoms_list)

        model = self._build_model()

        # NOTE: A real implementation would have a full training loop here.
        # For Cycle 1, we are just creating the model and saving it.
        # The actual training logic with mace-torch is complex and will
        # be expanded in a future cycle.

        torch.save(model, "trained_model.pt")

    def _create_z_table(self, atoms_list: list[Atoms]) -> AtomicNumberTable:
        """Creates a table of atomic numbers present in the dataset."""
        atomic_numbers = []
        for atoms in atoms_list:
            atomic_numbers.extend(atoms.get_atomic_numbers())
        return AtomicNumberTable([int(z) for z in sorted(list(set(atomic_numbers)))])

    def _prepare_data(self, atoms_list: list[Atoms]) -> list[AtomicData]:
        """Prepares the training data in the format required by MACE."""
        configurations = [config_from_atoms(atoms) for atoms in atoms_list]
        return [
            AtomicData.from_config(
                config, z_table=self.z_table, cutoff=self.config.r_cut
            )
            for config in configurations
        ]

    def _build_model(self) -> MACE:
        """Builds the MACE model based on the configuration."""
        # This is a simplified model configuration for demonstration purposes.
        model = MACE(
            r_max=self.config.r_cut,
            num_bessels=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls_first=None,
            interaction_cls=None,
            num_interactions=2,
            num_elements=len(self.z_table),
            hidden_irreps="128x0e + 128x1o",
            MLP_irreps="16x0e",
            atomic_energies=np.zeros(len(self.z_table)),  # Placeholder
            avg_num_neighbors=15,  # Placeholder
            atomic_numbers=self.z_table.zs,
            correlation=3,
        )
        return model
