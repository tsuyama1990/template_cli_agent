import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase.io import read, write

from mlip_autopipec.config import DFTInputConfig, DFTResult
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import ILabelingEngine, IProcessRunner


class LabelingEngine(ILabelingEngine):
    """
    Handles the execution of DFT calculations to label atomic structures.

    This class is responsible for taking an atomic structure from the database,
    generating the necessary input files for a Quantum Espresso calculation,
    running the calculation, parsing the output, and storing the results
    back in the database.
    """

    def __init__(
        self,
        dft_config: DFTInputConfig,
        db_wrapper: AseDBWrapper,
        process_runner: IProcessRunner,
        qe_command: str,
    ):
        """
        Initializes the LabelingEngine.

        Args:
            dft_config: Configuration for the DFT calculation.
            db_wrapper: Wrapper for the ASE database.
            process_runner: An object that implements the IProcessRunner interface.
            qe_command: The command to execute Quantum Espresso's pw.x.
        """
        self.dft_config = dft_config
        self.db_wrapper = db_wrapper
        self.process_runner = process_runner
        self.qe_command = qe_command

    def label_structure(self, structure_id: int) -> None:
        """
        Labels a single atomic structure by running a DFT calculation.

        Args:
            structure_id: The ID of the structure in the database to label.
        """
        atoms = self.db_wrapper.get_atoms_by_id(structure_id)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_file = temp_path / "qe_input.in"
            output_file = temp_path / "qe_output.out"

            write(
                input_file,
                atoms,
                format="espresso-in",
                pseudopotentials=self.dft_config.pseudopotentials,
                kpts=self.dft_config.kpoints,
                ecutwfc=self.dft_config.ecutwfc,
                control=self.dft_config.control,
            )

            try:
                command = [self.qe_command, "-in", str(input_file)]
                self.process_runner.run(command, str(output_file))

                result_atoms = read(output_file, format="espresso-out")
                energy = result_atoms.get_potential_energy()
                forces = result_atoms.get_forces()
                try:
                    stress = result_atoms.get_stress(voigt=False)
                except Exception:
                    stress = np.zeros((3, 3))

                dft_result = DFTResult(
                    energy=energy, forces=forces, stress=stress
                )
                self.db_wrapper.update_labels(structure_id, dft_result)

            except subprocess.CalledProcessError as e:
                print(f"Quantum Espresso execution failed: {e}")
                self.db_wrapper.update_state(
                    structure_id, "labeling_failed"
                )
