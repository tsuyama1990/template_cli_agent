import subprocess
from typing import Dict, Tuple

import ase
import numpy as np
from ase.io import write
from io import StringIO

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.utils.dft_utils import parse_qe_output


class LabellingEngine:
    """
    Encapsulates and manages the entire lifecycle of running a
    Quantum Espresso calculation.
    """

    def __init__(self, qe_command: str, pseudo_dir: str):
        self.qe_command = qe_command
        self.pseudo_dir = pseudo_dir

    def run(self, atoms: ase.Atoms) -> DFTResult:
        """
        Orchestrates the entire process for a single structure and returns
        a comprehensive DFTResult object.
        """
        try:
            input_content = self._generate_input_file(atoms)
            stdout, stderr, returncode = self._execute_qe(input_content)

            if returncode != 0 or "ERROR" in stdout or "ERROR" in stderr:
                return self._handle_error(atoms, stdout, stderr)

            parsed_data = parse_qe_output(stdout)

            if "energy" not in parsed_data or "forces" not in parsed_data:
                return self._handle_error(atoms, stdout, "Missing energy or forces")

            # Ensure stress is present, as it's required by the model
            if "stress" not in parsed_data:
                parsed_data["stress"] = np.zeros(6)

            return DFTResult(status="success", **parsed_data)

        except Exception as e:
            return self._handle_error(atoms, "", str(e))

    def _generate_input_file(self, atoms: ase.Atoms) -> str:
        """
        Translates an ase.Atoms object into a formatted Quantum Espresso input file.
        """
        # A simple example for a generic calculation.
        # This would be expanded in later cycles.
        pseudopotentials = {
            symbol: f"{symbol}.UPF" for symbol in set(atoms.get_chemical_symbols())
        }

        input_data = {
            "control": {
                "calculation": "scf",
                "restart_mode": "from_scratch",
                "prefix": "mlip",
                "outdir": "./out",
                "pseudo_dir": self.pseudo_dir,
            },
            "system": {
                "ibrav": 0,
                "ecutwfc": 30,
                "ecutrho": 240,
            },
            "electrons": {
                "conv_thr": 1e-8,
                "mixing_beta": 0.7,
            },
        }
        # In-memory write to a string
        f = StringIO()
        write(f, atoms, format="espresso-in", input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5)
        return f.getvalue()


    def _execute_qe(self, input_content: str) -> Tuple[str, str, int]:
        """
        Handles the execution of the pw.x binary in a secure subprocess.
        """
        process = subprocess.run(
            self.qe_command.split(),
            input=input_content,
            capture_output=True,
            text=True,
            check=False,
        )
        return process.stdout, process.stderr, process.returncode

    def _handle_error(self, atoms: ase.Atoms, stdout: str, stderr: str) -> DFTResult:
        """
        Implements the recovery logic for failed calculations.
        """
        return DFTResult(
            status="failed",
            energy=0.0,
            forces=np.zeros((len(atoms), 3)),
            stress=np.zeros(6),
            metadata={"stdout": stdout, "stderr": stderr},
        )
