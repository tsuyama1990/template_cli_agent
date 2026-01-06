"""
The core Molecular Dynamics (MD) and hybrid MD/MC simulation engine.
"""
from typing import List, Union
from ase import Atoms
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from ase.io import Trajectory
from ase.calculators.calculator import Calculator
from ase import units
import logging
from pathlib import Path

from mlip_autopipec.common.pydantic_models import MDConfig

# Configure logger for this module
logger = logging.getLogger(__name__)

class MDEngine:
    """
    Encapsulates the logic for running Molecular Dynamics simulations.

    For Cycle 1, this class provides basic NVT and NPT simulation capabilities.
    Advanced features like hybrid MD/MC and auto-ensemble switching are planned
    for future cycles.
    """

    def __init__(self, config: MDConfig, calculator: Calculator) -> None:
        """
        Initializes the MDEngine.

        Args:
            config: The MD configuration Pydantic model.
            calculator: An ASE-compatible calculator (e.g., MACE, EMT).
        """
        self.config = config
        self.calculator = calculator
        logger.info(f"Initialized MDEngine with config: {config}")

    def explore(self, initial_structures: List[Atoms], trajectory_path: Path) -> List[Atoms]:
        """
        Runs an MD simulation for each initial structure.

        Args:
            initial_structures: A list of starting ASE Atoms objects.
            trajectory_path: The path to save the output trajectory file.

        Returns:
            A list of all Atoms objects from the simulation trajectory.
        """
        full_trajectory = []
        for i, atoms in enumerate(initial_structures):
            logger.info(f"Starting MD exploration for structure {i+1}/{len(initial_structures)}...")

            atoms.calc = self.calculator

            dyn: Union[Langevin, NPT]
            # For Cycle 1, we use a simple logic: if pressure is near zero, use NVT.
            # Otherwise, use NPT. A more robust auto-detection is for Cycle 2.
            if abs(self.config.pressure_gpa) < 1e-5:
                dyn = self._setup_nvt(atoms)
                logger.info("Setting up NVT simulation (Langevin thermostat).")
            else:
                dyn = self._setup_npt(atoms)
                logger.info("Setting up NPT simulation.")

            # Attach a trajectory logger to the dynamics object
            with Trajectory(str(trajectory_path), 'a', atoms) as traj: # type: ignore
                dyn.attach(traj.write, interval=10) # type: ignore

                # Run the simulation for a fixed number of steps
                # This should be configurable in a real scenario
                dyn.run(200) # type: ignore

            # After the run, read the trajectory back to return it
            # In a real application, you might process this in a more memory-efficient way
            traj_atoms = Trajectory(str(trajectory_path)) # type: ignore
            full_trajectory.extend(list(traj_atoms))

            logger.info(f"MD exploration for structure {i+1} complete.")

        return full_trajectory

    def _setup_nvt(self, atoms: Atoms) -> Langevin:
        """Sets up an NVT simulation using the Langevin thermostat."""
        return Langevin(
            atoms,
            timestep=0.5 * units.fs,
            temperature_K=self.config.temperature_k,
            friction=0.02, # A reasonable default friction
        )

    def _setup_npt(self, atoms: Atoms) -> NPT:
        """Sets up an NPT simulation."""
        return NPT( # type: ignore
            atoms,
            timestep=0.5 * units.fs,
            temperature_K=self.config.temperature_k,
            externalstress=self.config.pressure_gpa * units.GPa,
            ttime=25 * units.fs, # Thermostat time constant
            pfactor=75**2 * units.fs**2 * units.GPa, # Barostat time constant
        )
