"""
This module implements the Explorer component of the MLIP-AutoPipe workflow.

The Explorer is responsible for taking a set of initial, static structures and
running high-speed, surrogate-based molecular dynamics (MD) simulations to
broadly sample the potential energy surface. It uses a pre-trained universal
potential, such as MACE, to make these simulations computationally inexpensive.
"""

import torch
from ase import Atoms
from ase.md.langevin import Langevin
from ase.units import fs
from mace.calculators import mace_mp

from ..config import ExplorerConfig, SimulationConfig
from ..database import AseDBWrapper
from ..interfaces import IExplorer


class Explorer(IExplorer):
    """
    Uses a surrogate model (MACE) to run MD simulations and explore phase space.
    """

    def __init__(
        self,
        explorer_config: ExplorerConfig,
        simulation_config: SimulationConfig,
        db_wrapper: AseDBWrapper,
    ):
        """
        Initializes the Explorer.

        Args:
            explorer_config: Configuration for the exploration phase.
            simulation_config: Configuration for simulation parameters (e.g., temperatures).
            db_wrapper: A wrapper for the ASE database to access initial structures.
        """
        self.config = explorer_config
        self.simulation_config = simulation_config
        self.db_wrapper = db_wrapper
        self.device = "cuda" if torch.cuda.is_available() else "cpu"

    def explore(self) -> list[Atoms]:
        """
        Runs the main exploration workflow.

        It retrieves initial structures, runs MD simulations for each one at
        various temperatures, and collects the resulting trajectories.

        Returns:
            A list of all atomic configurations sampled during the MD simulations.
        """
        initial_structures = self.db_wrapper.get_all_initial_structures()
        if not initial_structures:
            print("No initial structures found to explore.")
            return []

        print(
            f"Starting exploration with {len(initial_structures)} initial structures "
            f"at {len(self.simulation_config.temperature_steps)} temperatures."
        )

        all_explored_frames = []
        for i, atoms in enumerate(initial_structures):
            for temp in self.simulation_config.temperature_steps:
                print(f"Running MD for structure {i+1} at {temp} K...")
                frames = self._run_mace_md(atoms, temp)
                all_explored_frames.extend(frames)

        print(f"Exploration complete. Sampled {len(all_explored_frames)} total frames.")
        return all_explored_frames

    def _run_mace_md(self, atoms: Atoms, temperature: float) -> list[Atoms]:
        """
        Runs a single MD simulation for a given structure and temperature using MACE.

        Args:
            atoms: The initial atomic structure.
            temperature: The simulation temperature in Kelvin.

        Returns:
            A list of Atoms objects representing the MD trajectory.
        """
        if self.config.surrogate_model != "mace_mp":
            raise ValueError("Currently, only the 'mace_mp' surrogate model is supported.")

        # Load the pre-trained MACE-MP-0 model.
        # The model is automatically downloaded and cached by the MACE library.
        mace_calculator = mace_mp(
            model=None,  # Use default foundation model
            device=self.device,
            default_dtype="float64",
        )
        atoms.calc = mace_calculator

        # Set up the MD simulation using ASE's Langevin thermostat
        dyn = Langevin(
            atoms,
            timestep=2 * fs,
            temperature_K=temperature,
            friction=0.02,
        )

        # Attach an observer to collect the trajectory frames
        trajectory_frames = []

        def collect_frames():
            # We need to copy the atoms object to store the frame
            trajectory_frames.append(atoms.copy())

        # The interval determines how frequently we save a frame.
        # A value of 100 means we save a snapshot every 100 timesteps.
        dyn.attach(collect_frames, interval=100)

        # Run the simulation
        dyn.run(self.config.md_steps_per_structure)

        return trajectory_frames
