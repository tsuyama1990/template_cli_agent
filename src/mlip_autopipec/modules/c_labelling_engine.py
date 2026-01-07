import os
import subprocess
import tempfile

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils import dft_utils


class LabellingEngine:
    """
    Manages the process of labelling atomic structures using Quantum Espresso.
    """
    def __init__(self, qe_command: str, parameters: dict, db: AseDB):
        """
        Initializes the LabellingEngine.
        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "mpirun -np 4 pw.x").
            parameters: A dictionary of parameters for the QE calculation.
            db: An instance of the AseDB wrapper.
        """
        if not qe_command:
            raise ValueError("Quantum Espresso command cannot be empty.")
        self._qe_command = qe_command.split()
        self._parameters = parameters
        self._db = db

    def execute(self, atoms: Atoms) -> int:
        """
        Generates input, runs a QE calculation, parses the output, and saves the result.
        The process is executed in a temporary directory to keep files isolated.

        Args:
            atoms: The ASE Atoms object to be calculated.

        Returns:
            The database ID of the new entry.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            input_content = dft_utils.generate_qe_input(atoms, self._parameters)

            input_filepath = os.path.join(temp_dir, "qe.in")
            output_filepath = os.path.join(temp_dir, "qe.out")

            with open(input_filepath, "w") as f:
                f.write(input_content)

            # Construct the full command
            command = self._qe_command + ["-in", input_filepath]

            # Execute the command and write stdout to the output file
            process_result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                cwd=temp_dir
            )

            with open(output_filepath, "w") as f:
                f.write(process_result.stdout)
                if process_result.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(process_result.stderr)

            # Parse the output
            dft_result = dft_utils.parse_qe_output(process_result.stdout)

            # Write to database
            db_id = self._db.write(atoms, dft_result)

        return db_id
