import os
import tempfile

from ase import Atoms
from ase.calculators.emt import EMT
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.units import fs

from mlip_autopipec.common.pydantic_models import FullConfig


class MDEngine:
    """
    Runs Molecular Dynamics simulations using ASE, designed for robustness
    and compatibility with multiprocessing.
    """

    def __init__(self, config: FullConfig):
        self.config = config.exploration

    def _get_calculator(self):
        """
        Instantiates and returns the ASE calculator.
        This 'Late Binding' is crucial for multiprocessing, as it prevents
        pickling large, complex calculator objects (like MLIP models).
        """
        # For Cycle 1, we use a simple, fast calculator. In a real scenario,
        # this would load a heavyweight MLIP model (e.g., MACE, SevenNet).
        return EMT()

    def run(self, atoms: Atoms) -> list[Atoms]:
        """
        Executes an MD simulation for a given initial structure.

        Args:
            atoms: An ASE Atoms object representing the initial structure.

        Returns:
            A list of Atoms objects representing the MD trajectory.
        """
        # 1. Late Binding: Instantiate and attach the calculator inside the run method.
        atoms.calc = self._get_calculator()

        # Use a temporary file to store the trajectory, ensuring data is written to disk
        # before being processed by the next pipeline stage.
        with tempfile.NamedTemporaryFile(suffix=".traj", delete=False) as tmp:
            traj_path = tmp.name

        try:
            # 2. Set up the ASE dynamics engine (Langevin for NVT ensemble)
            dyn = Langevin(
                atoms,
                timestep=self.config.time_step_fs * fs,
                temperature_K=self.config.temperature_k,
                friction=0.02,  # A reasonable default friction value
            )

            # Attach a logger to save the trajectory at specified intervals
            traj = Trajectory(traj_path, 'w', atoms)
            # Save every 10 steps for efficiency in this example
            dyn.attach(traj.write, interval=max(1, self.config.total_steps // 10))

            # 3. Run the simulation.
            # Robust implementation would wrap this in a try/except block to
            # catch ASE-specific errors from simulations that fail, allowing
            # the pipeline to continue with other structures.
            dyn.run(self.config.total_steps)

            # 4. Read the full trajectory back from the file.
            traj_out = Trajectory(traj_path)
            return [at for at in traj_out]

        finally:
            # 5. Ensure cleanup of the temporary file in all cases.
            if os.path.exists(traj_path):
                os.remove(traj_path)
