"""The Labelling Engine for performing automated DFT calculations."""
import subprocess
import tempfile
import os
from typing import Dict, Any
from ase.atoms import Atoms
from ..data.database import AseDB
from ..utils.dft_utils import generate_qe_input, parse_qe_output

class LabellingEngine:
    """Orchestrates the DFT calculation for a given atomic structure."""

    def __init__(self, qe_command: str, parameters: Dict[str, Any], pseudos: Dict[str, str], db: AseDB):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso's pw.x (e.g., "mpirun -np 4 pw.x").
            parameters: Dictionary of QE input parameters.
            pseudos: Dictionary mapping atom species to pseudopotential filenames.
            db: An instance of the AseDB wrapper.
        """
        self._qe_command = qe_command
        self._parameters = parameters
        self._pseudos = pseudos
        self._db = db

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation for the given Atoms object and stores the result in the database.

        Args:
            atoms: The ASE Atoms object to be calculated.

        Returns:
            The ID of the new entry in the database.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "qe.in")
            output_path = os.path.join(temp_dir, "qe.out")

            # Generate input file
            input_content = generate_qe_input(atoms, self._parameters, self._pseudos)
            with open(input_path, "w") as f:
                f.write(input_content)

            # Construct and run the command
            command = f"{self._qe_command} -in {input_path}"

            completed_process = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                cwd=temp_dir
            )

            # Parse the output
            output_content = completed_process.stdout
            if completed_process.returncode != 0:
                # Append stderr for more context on failure
                output_content += "\n--- STDERR ---\n" + completed_process.stderr

            dft_result = parse_qe_output(output_content)

            # Write to database
            db_id = self._db.write(atoms, dft_result)
            return db_id
