"""Application services for the MLIP AutoPipe workflow."""

import io
import torch
import numpy as np
from ase import Atoms
from ase.io import write
from ase.calculators.lj import LennardJones
from mace.data.atomic_data import AtomicData, get_data_loader
from mace.data.utils import config_from_atoms
from mace.tools import AtomicNumberTable
from mace.modules.blocks import RealAgnosticInteractionBlock, InteractionBlock
from mace.modules.loss import WeightedEnergyForcesLoss
from mace.modules.models import MACE
from e3nn.o3 import Irreps

from mlip_autopipec.domain.ports import DatabasePort, ProcessRunnerPort
from mlip_autopipec.infrastructure import dft_utils
from mlip_autopipec.config import Settings

class LabellingService:
    """Orchestrates the DFT labelling of atomic structures."""
    def __init__(self, db: DatabasePort, runner: ProcessRunnerPort, settings: Settings):
        self.db = db
        self.runner = runner
        self.settings = settings

    def run(self):
        unlabelled_atoms = self.db.get_atoms_by_state('unlabelled')
        print(f"Found {len(unlabelled_atoms)} unlabelled structures.")
        for atoms in unlabelled_atoms:
            try:
                input_str = self._generate_qe_input(atoms)
                output_str = self._execute_dft(input_str)
                dft_result = dft_utils.parse_qe_output(output_str)
                self.db.update_with_dft_results(atoms.info['db_id'], dft_result)
            except Exception as e:
                print(f"Error processing structure ID {atoms.info['db_id']}: {e}")

    def _generate_qe_input(self, atoms: Atoms) -> str:
        input_data = {
            'control': {
                'calculation': 'scf', 'prefix': 'pwscf',
                'pseudo_dir': self.settings.pseudo_dir, 'outdir': './out',
            },
            'system': {
                'ibrav': 0, 'nat': len(atoms),
                'ntyp': len(set(atoms.get_chemical_symbols())),
                'ecutwfc': self.settings.ecutwfc,
            },
            'electrons': {'mixing_beta': 0.7},
        }
        string_io = io.StringIO()
        write(
            string_io, atoms, format='espresso-in', input_data=input_data,
            pseudopotentials={s: 'pseudo.UPF' for s in set(atoms.get_chemical_symbols())},
            kpts=(1, 1, 1),
        )
        return string_io.getvalue()

    def _execute_dft(self, input_str: str) -> str:
        command = self.settings.dft_command.split()
        return self.runner.run(command, input_str)

class TrainingService:
    """Orchestrates the training of the MLIP model."""
    def __init__(self, db: DatabasePort, settings: Settings):
        self.db = db
        self.settings = settings
        self.z_table = None

    def run(self):
        configs, dft_results = self._prepare_dataset()
        if not configs:
            return
        baseline_results = self._calculate_baseline(configs)
        delta_results = self._calculate_delta(dft_results, baseline_results)
        self._train_model(configs, delta_results)

    def _prepare_dataset(self) -> (list, list):
        labelled_atoms = self.db.get_atoms_by_state('labelled')
        atomic_numbers = {n for a in labelled_atoms for n in a.get_atomic_numbers()}
        if not atomic_numbers:
            return [], []
        self.z_table = AtomicNumberTable([int(z) for z in sorted(list(atomic_numbers))])
        configs = [config_from_atoms(atoms) for atoms in labelled_atoms]
        dft_results = [{'energy': a.get_potential_energy(), 'forces': a.get_forces()} for a in labelled_atoms]
        return configs, dft_results

    def _calculate_baseline(self, configs: list) -> list:
        lj = LennardJones(sigma=3.4, epsilon=0.0103)
        baseline_results = []
        for config in configs:
            atoms = Atoms(
                numbers=config.atomic_numbers, positions=config.positions,
                cell=config.cell, pbc=config.pbc,
            )
            atoms.calc = lj
            baseline_results.append({'energy': atoms.get_potential_energy(), 'forces': atoms.get_forces()})
        return baseline_results

    def _calculate_delta(self, dft_results: list, baseline_results: list) -> list:
        delta_results = []
        for dft, baseline in zip(dft_results, baseline_results, strict=True):
            delta_results.append({
                'energy': dft['energy'] - baseline['energy'],
                'forces': dft['forces'] - baseline['forces']
            })
        return delta_results

    def _train_model(self, configs: list, delta_results: list):
        atomic_data = [
            AtomicData.from_config(c, z_table=self.z_table, cutoff=self.settings.cutoff) for c in configs
        ]
        for ad, delta in zip(atomic_data, delta_results, strict=True):
            ad.energy = torch.tensor(delta['energy'], dtype=torch.float32)
            ad.forces = torch.tensor(delta['forces'], dtype=torch.float32)

        model_config = {
            'r_max': self.settings.cutoff, 'num_bessel': 8, 'num_polynomial_cutoff': 5,
            'max_ell': 3, 'interaction_cls': RealAgnosticInteractionBlock,
            'interaction_cls_first': InteractionBlock,
            'num_interactions': 2, 'hidden_irreps': Irreps('128x0e + 128x1o'),
            'MLP_irreps': Irreps('16x0e'),
            'gate': torch.nn.functional.silu,
            'atomic_energies': np.zeros(len(self.z_table)),
            'avg_num_neighbors': 15.0,
            'atomic_numbers': self.z_table.zs,
            'correlation': 3,
            'num_elements': len(self.z_table.zs)
        }
        model = MACE(**model_config)

        loss_fn = WeightedEnergyForcesLoss(energy_weight=1.0, forces_weight=1.0)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

        data_loader = get_data_loader(atomic_data, batch_size=1)

        print("Starting training loop...")
        for epoch in range(2):
            for batch in data_loader:
                optimizer.zero_grad()
                output = model(batch)
                loss = loss_fn(output, batch)
                loss.backward()
                optimizer.step()
            print(f"Epoch {epoch+1}, Loss: {loss.item()}")

        torch.save(model.state_dict(), self.settings.model_path)
        print(f"Model saved to {self.settings.model_path}")
