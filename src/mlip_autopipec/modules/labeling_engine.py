import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from mlip_autopipec.config import DFTInputConfig, DFTResult
from mlip_autopipec.interfaces import ILabelingEngine, IProcessRunner


class LabelingEngine(ILabelingEngine):
    """
    Handles the execution of DFT calculations to label atomic structures.

    This class is responsible for taking an atomic structure, generating the
    necessary input files for a Quantum Espresso calculation, running the
    calculation, parsing the output, and returning the results.
    """

    def __init__(
        self,
        dft_input_configuration: DFTInputConfig,
        process_runner: IProcessRunner,
        qe_command: str,
    ):
        """
        Initializes the LabelingEngine.

        Args:
            dft_input_configuration: Configuration for the DFT calculation.
            process_runner: An object that implements the IProcessRunner interface.
            qe_command: The command to execute Quantum Espresso's pw.x.
        """
        self.dft_input_configuration = dft_input_configuration
        self.process_runner = process_runner
        self.qe_command = qe_command

    def label_structure(self, atoms: Atoms) -> DFTResult:
        """
        Labels a single atomic structure by running a DFT calculation.

        Args:
            atoms: The ASE Atoms object to label.

        Returns:
            A DFTResult object containing the calculated energy, forces, and stress.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_file = temp_path / "qe_input.in"
            output_file = temp_path / "qe_output.out"

            write(
                input_file,
                atoms,
                format="espresso-in",
                pseudopotentials=self.dft_input_configuration.pseudopotentials,
                kpts=self.dft_input_configuration.kpoints,
                ecutwfc=self.dft_input_configuration.ecutwfc,
                control=self.dft_input_configuration.control,
            )

            command = [self.qe_command, "-in", str(input_file)]
            self.process_runner.run(command, str(output_file))

            result_atoms = read(output_file, format="espresso-out")
            energy = result_atoms.get_potential_energy()
            forces = result_atoms.get_forces()
            try:
                stress = result_atoms.get_stress(voigt=False)
            except Exception:
                stress = np.zeros((3, 3))

            return DFTResult(energy=energy, forces=forces, stress=stress)
