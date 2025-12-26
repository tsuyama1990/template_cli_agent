import subprocess
import tempfile
import os
from pathlib import Path

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.utils import dft_utils

class LabellingEngine:
    """
    Manages the execution of DFT calculations using Quantum Espresso.
    """

    def __init__(self, qe_command: str, db: AseDB):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x").
            db: An instance of the AseDB class for database operations.
        """
        self._qe_command = qe_command
        self._db = db

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a Quantum Espresso calculation for a given atomic structure.

        Generates input, runs QE in a temporary directory, parses the output,
        and saves the result to the database.

        Args:
            atoms: The atomic structure to calculate.

        Returns:
            The database ID of the new entry.
        """
        result: DFTResult
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            input_path = tmpdir_path / "qe.in"
            output_path = tmpdir_path / "qe.out"

            # 1. Generate QE input file
            qe_input = dft_utils.generate_qe_input(atoms)
            with open(input_path, "w") as f:
                f.write(qe_input)

            # 2. Execute Quantum Espresso
            command = f"{self._qe_command} -in {input_path}"
            process = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                cwd=tmpdir_path,
            )

            # 3. Parse the output
            if process.returncode == 0:
                result = dft_utils.parse_qe_output(process.stdout)
            else:
                error_msg = f"Quantum Espresso exited with code {process.returncode}.\nStderr: {process.stderr}"
                result = DFTResult(
                    total_energy_ev=0.0,
                    forces=[],
                    stress=[],
                    was_successful=False,
                    error_message=error_msg,
                )

        # 4. Save to database
        db_id = self._db.write(atoms, result)
        return db_id
