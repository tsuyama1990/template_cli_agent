import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import kpts2mp
from ase.io import read, write

from mlip_autopipec.config import DFTComputeConfig, DFTResult
from mlip_autopipec.interfaces import ILabelingEngine, IProcessRunner


class LabelingEngine(ILabelingEngine):
    """
    Handles the execution of DFT calculations to label atomic structures.
    """

    def __init__(
        self,
        dft_compute_config: DFTComputeConfig,
        process_runner: IProcessRunner,
        qe_command: str,
    ):
        """Initializes the LabelingEngine."""
        self.dft_compute_config = dft_compute_config
        self.process_runner = process_runner
        self.qe_command = qe_command

    def label_structure(self, atoms: Atoms) -> DFTResult:
        """Labels a single atomic structure by running a DFT calculation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_file = temp_path / "qe_input.in"
            output_file = temp_path / "qe_output.out"

            kpts_array = kpts2mp(atoms, self.dft_compute_config.kpoints_density, even=True)
            kpts = tuple(kpts_array)

            write(
                input_file,
                atoms,
                format="espresso-in",
                pseudopotentials=self.dft_compute_config.pseudopotentials,
                ecutwfc=self.dft_compute_config.ecutwfc,
                ecutrho=self.dft_compute_config.ecutrho,
                kpts=kpts,
                control=self.dft_compute_config.control,
            )

            command = [self.qe_command, "-in", str(input_file)]
            self.process_runner.run(command, str(output_file))

            result_atoms = read(output_file, format="espresso-out")

            # Ensure calculator is attached for property access
            if not result_atoms.calc:
                raise RuntimeError("ASE failed to parse calculator results from QE output.")

            energy = result_atoms.get_potential_energy()
            forces = result_atoms.get_forces()
            try:
                stress = result_atoms.get_stress(voigt=False)
            except (IndexError, ValueError):  # Handle ASE parsing errors
                stress = np.zeros((3, 3))

            return DFTResult(energy=energy, forces=forces, stress=stress)
