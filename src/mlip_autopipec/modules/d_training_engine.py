import logging

import numpy as np
import torch
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.calculators.singlepoint import SinglePointCalculator
from mace.data.atomic_data import AtomicData
from mace.data.utils import config_from_atoms
from mace.modules.models import MACE
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import MLIPTraining

Z_TABLE = AtomicNumberTable([1, 6, 8, 14])
logger = logging.getLogger(__name__)


class TrainingEngine:
    """
    Manages the process of training an MLIP model on labeled data.
    """

    def __init__(self, config: MLIPTraining, db_wrapper: AseDBWrapper):
        """
        Initializes the TrainingEngine.

        Args:
            config: The MLIPTraining configuration object.
            db_wrapper: An instance of the AseDBWrapper.
        """
        self.config = config
        self.db_wrapper = db_wrapper

    def execute(self):
        """
        Executes the training workflow.
        """
        logger.info("Starting Training Engine...")
        labeled_rows = self.db_wrapper.get_all_labeled_rows()
        if not labeled_rows:
            logger.info("No labeled structures found to train on.")
            return

        logger.info(f"Found {len(labeled_rows)} labeled structures for training.")
        atoms_list = [row.toatoms() for row in labeled_rows]
        training_data = self._prepare_training_data(atoms_list)

        model = MACE(
            r_max=self.config.r_cut,
            num_bessel=8,
            num_polynomial_basis=4,
            radial_MLP=[64, 64, 64],
            interaction_cls_first=None,
            interaction_cls=None,
            hidden_irreps="128x0e + 128x1o",
            atomic_energies=np.zeros(len(Z_TABLE.zs)),
            avg_num_neighbors=15,
            atomic_numbers=Z_TABLE.zs,
            correlation=3,
            gate=torch.nn.functional.silu,
        )

        logger.info("Starting simplified training loop...")
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        loss_fn = torch.nn.MSELoss()

        for epoch in range(5):
            for data in training_data:
                optimizer.zero_grad()
                config = config_from_atoms(data)
                atomic_data = AtomicData.from_config(
                    config, z_table=Z_TABLE, cutoff=self.config.r_cut
                )
                output = model(atomic_data.to_dict())
                target_energy = torch.tensor(
                    [data.get_potential_energy()], dtype=torch.float32
                )
                loss = loss_fn(output['energy'], target_energy)
                loss.backward()
                optimizer.step()
            logger.info(f"  Epoch {epoch+1}/5, Loss: {loss.item():.6f}")

        model_path = "trained_model.pt"
        torch.save(model.state_dict(), model_path)
        logger.info(f"Training finished. Model saved to '{model_path}'.")

    def _prepare_training_data(self, atoms_list: list[Atoms]) -> list[Atoms]:
        """
        Prepares the training data, applying delta learning if enabled.
        """
        if not self.config.delta_learning:
            return atoms_list

        logger.info(
            f"Applying delta learning with base potential: {self.config.base_potential}"
        )
        if self.config.base_potential != 'lj_auto':
            raise NotImplementedError(
                "Only 'lj_auto' is supported for base_potential in Cycle 1."
            )

        base_calc = LennardJones()
        processed_atoms_list = []

        for atoms in atoms_list:
            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            temp_atoms = atoms.copy()
            temp_atoms.calc = base_calc
            base_energy = temp_atoms.get_potential_energy()
            base_forces = temp_atoms.get_forces()

            delta_energy = dft_energy - base_energy
            delta_forces = dft_forces - base_forces

            delta_atoms = atoms.copy()
            delta_atoms.calc = SinglePointCalculator(
                delta_atoms, energy=delta_energy, forces=delta_forces
            )

            processed_atoms_list.append(delta_atoms)

        return processed_atoms_list
