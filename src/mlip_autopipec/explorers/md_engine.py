# src/mlip_autopipec/explorers/md_engine.py
"""
Encapsulates the logic for running Molecular Dynamics (MD) simulations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms, units
from ase.calculators.calculator import Calculator
from ase.md.langevin import Langevin

if TYPE_CHECKING:
    from mlip_autopipec.common.pydantic_models import FullConfig


class MDEngine:
    """Runs a Molecular Dynamics simulation using ASE."""

    def __init__(self, config: FullConfig, calculator: Calculator) -> None:
        """
        Initializes the MD engine.

        Args:
            config: The full application configuration.
            calculator: An ASE-compatible calculator (e.g., EMT, MACE).
        """
        self.config = config
        self.calculator = calculator

    def explore(self, atoms: Atoms) -> list[Atoms]:
        """
        Runs the MD simulation on a given initial structure.

        Args:
            atoms: The starting atomic structure.

        Returns:
            A list of Atoms objects representing the MD trajectory.
        """
        atoms.set_calculator(self.calculator)

        md_config = self.config.exploration

        # Use Langevin dynamics for NVT ensemble
        dyn = Langevin(
            atoms,
            md_config.timestep_fs * units.fs,
            temperature_K=md_config.temperature_k,
            friction=0.02,
        )

        trajectory = []

        # Callback function to store each frame of the trajectory
        def log_step() -> None:
            # Create a copy to store the state at this step
            frame = atoms.copy()  # type: ignore[no-untyped-call]
            # The calculator is not automatically copied, so re-attach it
            frame.calc = atoms.calc
            frame.info["energy"] = frame.get_potential_energy()
            frame.arrays["forces"] = frame.get_forces()
            trajectory.append(frame)

        dyn.attach(log_step, interval=1)  # type: ignore[no-untyped-call]

        # Run the dynamics
        dyn.run(md_config.n_steps)  # type: ignore[no-untyped-call]

        return trajectory
