from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from ase.calculators.emt import EMT
from ase.md.langevin import Langevin

# from ase.md.mc.swap import SwapMC
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import LBFGS

from mlip_autopipec.core.physics import detect_vacuum

if TYPE_CHECKING:
    from ase import Atoms

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MDExplorer:
    """
    Computational engine for exploring the potential energy surface.

    This class takes a set of seed structures and uses Molecular Dynamics (MD)
    and hybrid Monte Carlo (MC) simulations to generate a diverse ensemble of
    atomic configurations.
    """

    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.md_steps = config.get("md_steps", 1000)
        self.temperature = config.get("temperature", 300.0)
        self.pressure = config.get("pressure", 1.0)  # In GPa
        self.mlip_model = config.get("mlip_model", "emt")
        self.use_zbl = config.get("use_zbl", True)
        self.swap_frequency = config.get("swap_frequency", 10)

    def _get_calculator(self):
        """Lazily loads the calculator to improve performance."""
        if self.mlip_model == "emt":
            mlip_calc = EMT()
        else:
            # In a real scenario, this would load a MACE or NequIP model
            # from ase.calculators.mace import MACECalculator
            # calculator = MACECalculator(model_path=self.mlip_model)
            raise NotImplementedError(f"MLIP model '{self.mlip_model}' is not supported yet.")

        return mlip_calc

    def explore(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Run MD/MC simulations on a list of structures.

        Args:
            structures: A list of seed `ase.Atoms` objects.

        Returns:
            A list of new `ase.Atoms` objects from the simulation trajectories.
        """
        all_frames = []
        for i, atoms in enumerate(structures):
            logger.info(f"Exploring structure {i + 1}/{len(structures)}...")
            try:
                atoms.calc = self._get_calculator()
                # Initial relaxation before MD
                opt = LBFGS(atoms, trajectory=None, logfile=None)
                opt.run(fmax=0.1, steps=50)

                MaxwellBoltzmannDistribution(atoms, temperature_K=self.temperature)

                ensemble_type = detect_vacuum(atoms)
                logger.info(f"Detected ensemble type: {ensemble_type}")

                if ensemble_type == "bulk":
                    dyn = NPT(
                        atoms,
                        timestep=1.0,
                        temperature_K=self.temperature,
                        external_pressure=self.pressure,
                        ttime=25.0,
                        pfactor=75.0,
                    )
                else:  # Slab or molecule
                    dyn = Langevin(
                        atoms,
                        timestep=1.0,
                        temperature_K=self.temperature,
                        friction=0.02,
                    )


                def collect_frames():
                    all_frames.append(atoms.copy())

                dyn.attach(collect_frames, interval=self.md_steps // 10)
                dyn.run(self.md_steps)

            except Exception as e:
                logger.error(
                    f"MD simulation failed for structure {i}: {e}",
                    exc_info=True,
                )
                continue  # Move to the next structure

        return all_frames
