import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase.io import read, write

from mlip_autopipec.config import DFTInputConfig, DFTResult
from mlip_autopipec.database import AseDBWrapper


class LabelingEngine:
    """Handles the execution of DFT calculations to label atomic structures."""

    def __init__(
        self,
        dft_config: DFTInputConfig,
        db_wrapper: AseDBWrapper,
        qe_command: str,
    ):
        """Initializes the LabelingEngine.

        Args:
            dft_config: Configuration for the DFT calculation.
            db_wrapper: Wrapper for the ASE database.
            qe_command: The command to execute Quantum Espresso's pw.x.
        """
        self.dft_config = dft_config
        self.db_wrapper = db_wrapper
        self.qe_command = qe_command

    def label_structure(self, id: int):
        """Labels a single atomic structure by running a DFT calculation.

        Args:
            id: The ID of the structure in the database to label.
        """
        atoms = self.db_wrapper.get_atoms_by_id(id)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_file = temp_path / "qe_input.in"
            output_file = temp_path / "qe_output.out"

            # Write QE input file
            write(
                input_file,
                atoms,
                format="espresso-in",
                pseudopotentials=self.dft_config.pseudopotentials,
                kpts=self.dft_config.kpoints,
                ecutwfc=self.dft_config.ecutwfc,
                control=self.dft_config.control,
            )

            # Execute QE
            command = [self.qe_command, "-in", str(input_file)]
            with open(output_file, "w") as f_out:
                subprocess.run(  # noqa: S603
                    command, stdout=f_out, shell=False, check=True
                )

            # Parse QE output
            result_atoms = read(output_file, format="espresso-out")
            energy = result_atoms.get_potential_energy()
            forces = result_atoms.get_forces()
            try:
                stress = result_atoms.get_stress(voigt=False)
            except Exception:
                # Stress parsing can fail on some QE versions/setups, default to zero.
                stress = np.zeros((3, 3))

            dft_result = DFTResult(energy=energy, forces=forces, stress=stress)
            self.db_wrapper.update_labels(id, dft_result)
