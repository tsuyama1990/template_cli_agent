
from ase import Atoms
from ase.calculators.emt import EMT
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet

from mlip_autopipec.common.pydantic_models import ExplorationConfig, MDEnsemble
from mlip_autopipec.interfaces import IExplorer


class MDEngine(IExplorer):
    def __init__(self, config: ExplorationConfig):
        self.config = config

    def run(self, initial_structures: list[Atoms]) -> list[Atoms]:
        """
        Runs an MD simulation for each initial structure and returns the trajectory.
        """
        all_frames = []
        for atoms in initial_structures:
            # Assign a calculator - using EMT as a placeholder for a real MLIP
            # In a real scenario, this would load the specified calculator (e.g., MACE)
            atoms.calc = EMT()

            # Set initial temperature
            MaxwellBoltzmannDistribution(atoms, temperature_K=self.config.temperature_k)

            # Select dynamics ensemble
            if self.config.ensemble == MDEnsemble.nvt:
                dyn = VelocityVerlet(atoms, timestep=self.config.time_step_fs)
            elif self.config.ensemble == MDEnsemble.npt:
                dyn = NPT(
                    atoms,
                    timestep=self.config.time_step_fs,
                    temperature_K=self.config.temperature_k,
                    externalstress=0.0,
                )
            else:
                raise ValueError(f"Unknown ensemble: {self.config.ensemble}")

            def collect_frames():
                all_frames.append(atoms.copy())

            dyn.attach(collect_frames, interval=1)
            dyn.run(self.config.num_steps)

        return all_frames
