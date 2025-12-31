from typing import Any, List, Tuple

import ase
import numpy as np
import torch
import mace
from ase.calculators.lj import LennardJones
from mace.data.atomic_data import AtomicData
from mace.data.utils import config_from_atoms
from mace.tools import AtomicNumberTable

from mlip_autopipec.data.models import DFTResult


class TrainingEngine:
    """
    Abstracts the complexities of the underlying MLIP framework and
    implements the delta learning logic.
    """

    def __init__(self, model_type: str = "mace"):
        if model_type != "mace":
            raise NotImplementedError("Only MACE is supported in Cycle 1.")
        self.model_type = model_type
        self.z_table = None

    def train(self, data: List[Tuple[ase.Atoms, DFTResult]]) -> Any:
        """
        Takes the full labelled dataset and returns a trained,
        serialisable model object.
        """
        if not data:
            raise ValueError("Training data cannot be empty.")

        # 1. Prepare dataset for MACE
        all_atoms = [d[0] for d in data]
        all_elements = set()
        for atoms in all_atoms:
            all_elements.update(atoms.get_chemical_symbols())

        self.z_table = AtomicNumberTable([ase.data.atomic_numbers[s] for s in sorted(list(all_elements))])

        baseline_potential = self._get_baseline_potential(list(all_elements))

        mace_dataset = []
        for atoms, dft_result in data:
            # 2. Compute delta
            delta_energy, delta_forces = self._compute_delta(
                atoms, dft_result, baseline_potential
            )

            # 3. Create MACE-compatible data object
            config = config_from_atoms(atoms)
            atomic_data = AtomicData.from_config(
                config, z_table=self.z_table, cutoff=3.0 # A default cutoff
            )
            atomic_data.energy = torch.tensor(delta_energy).unsqueeze(0)
            atomic_data.forces = torch.tensor(delta_forces)
            mace_dataset.append(atomic_data)

        # 4. Simple MACE model for demonstration
        # In a real scenario, this would be a more complex setup
        from mace.modules.models import MACE
        model = MACE(
            r_max=3.0,
            num_bessel=8,
            num_polynomial_cutoff=5,
            max_ell=3,
            interaction_cls=mace.modules.blocks.RealAgnosticInteractionBlock,
            interaction_cls_first=mace.modules.blocks.RealAgnosticInteractionBlock,
            num_interactions=2,
            num_elements=len(all_elements),
            hidden_irreps="16x0e",
            atomic_energies=torch.zeros(len(all_elements)), # Not used in delta learning
            avg_num_neighbors=8, # A reasonable default
            atomic_numbers=self.z_table.zs
        )

        # 5. Dummy training loop
        # A real implementation would use a proper optimizer and loss function
        # and iterate over the data multiple times.
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
        for atomic_data in mace_dataset:
            optimizer.zero_grad()
            output = model(atomic_data.to_dict())
            loss = (output["energy"] - atomic_data.energy).pow(2).sum()
            loss += (output["forces"] - atomic_data.forces).pow(2).sum()
            loss.backward()
            optimizer.step()

        return model


    def _get_baseline_potential(
        self, elements: List[str]
    ) -> ase.calculators.calculator.Calculator:
        """
        Returns a simple, physics-based ASE calculator for delta learning.
        """
        # Using default LJ parameters from ASE, which is fine for a baseline.
        return LennardJones()


    def _compute_delta(
        self,
        atoms: ase.Atoms,
        dft_result: DFTResult,
        baseline: ase.calculators.calculator.Calculator,
    ) -> Tuple[float, np.ndarray]:
        """
        Calculates the residual between the DFT ground truth and the baseline potential.
        """
        atoms.calc = baseline
        baseline_energy = atoms.get_potential_energy()
        baseline_forces = atoms.get_forces()

        delta_energy = dft_result.energy - baseline_energy
        delta_forces = dft_result.forces - baseline_forces

        return delta_energy, delta_forces
