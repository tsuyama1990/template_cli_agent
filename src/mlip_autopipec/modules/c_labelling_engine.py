import os
import subprocess
import tempfile

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output


class LabellingEngine:
    """
    The Labelling Engine (Module C) for performing automated DFT calculations.
    """
    def __init__(self, qe_command: str, db: AseDB, default_params: dict = None):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x").
            db: An instance of the AseDB wrapper.
            default_params: Default parameters for QE calculations.
        """
        self._qe_command = qe_command
        self._db = db
        self._default_params = default_params or {}

    def execute(self, atoms: Atoms) -> int:
        """
        Takes an ASE Atoms object, runs a DFT calculation, and stores the result.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            qe_input = generate_qe_input(atoms, self._default_params)
            input_path = os.path.join(tmpdir, "QE.in")
            output_path = os.path.join(tmpdir, "QE.out")

            with open(input_path, "w") as f:
                f.write(qe_input)

            command = [self._qe_command, "-in", input_path]

            try:
                with open(output_path, "w") as outfile:
                    subprocess.run(
                        command,
                        cwd=tmpdir,
                        check=True,
                        capture_output=True,
                        text=True,
                        stdout=outfile,
                        stderr=subprocess.PIPE,
                    )

                with open(output_path) as f:
                    qe_output = f.read()

                energy, forces, stress, success, error = parse_qe_output(qe_output)

                if not success:
                    dft_result = DFTResult(
                        total_energy_ev=0, forces=[], stress=[], was_successful=False,
                        error_message=error
                    )
                else:
                    calc = SinglePointCalculator(
                        atoms, energy=energy, forces=forces, stress=stress
                    )
                    atoms.calc = calc
                    dft_result = DFTResult(
                        total_energy_ev=energy, forces=forces, stress=stress,
                        was_successful=True
                    )

            except subprocess.CalledProcessError as e:
                error_msg = (
                    f"Quantum Espresso execution failed with exit code {e.returncode}.\n"
                    f"Stderr: {e.stderr}"
                )
                dft_result = DFTResult(
                    total_energy_ev=0, forces=[], stress=[], was_successful=False,
                    error_message=error_msg
                )
            except FileNotFoundError:
                error_msg = (
                    f"Command '{self._qe_command}' not found. "
                    "Is Quantum Espresso in your PATH?"
                )
                dft_result = DFTResult(
                    total_energy_ev=0, forces=[], stress=[], was_successful=False,
                    error_message=error_msg
                )

        db_id = self._db.write(atoms, dft_result.model_dump())
        return db_id
