# Description: The Labelling Engine (Module C) for automated DFT calculations.
import shutil
import subprocess
import uuid
from pathlib import Path

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output


class LabellingEngine:
    """Orchestrates the DFT calculation process using Quantum Espresso."""

    def __init__(self, qe_command: str, db: AseDB):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x -in QE.in").
            db: An instance of the AseDB wrapper.
        """
        self._qe_command = qe_command
        self._db = db

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation for a given atomic structure and saves the result.

        This method creates a temporary directory, generates the QE input file,
        runs the calculation, parses the output, and stores the result in the
        database before cleaning up the temporary directory.

        Args:
            atoms: The ase.Atoms object to be calculated.

        Returns:
            The database ID of the newly created entry.
        """
        # Create a unique temporary directory for the calculation
        temp_dir = Path(f"./temp_qe_run_{uuid.uuid4()}")
        temp_dir.mkdir(parents=True, exist_ok=True)
        out_dir = temp_dir / "out"
        out_dir.mkdir(parents=True, exist_ok=True)

        input_filename = "QE_input.in"
        output_filename = "QE_output.out"
        input_filepath = temp_dir / input_filename
        output_filepath = temp_dir / output_filename

        try:
            # Ensure the parent directory exists
            input_filepath.parent.mkdir(parents=True, exist_ok=True)

            # Generate and write the QE input file
            input_str = generate_qe_input(atoms, {"CONTROL": {"outdir": f"'{out_dir}/'"}})
            input_filepath.write_text(input_str)

            # Construct and run the command
            # Using a list of args is safer than shell=True
            command_parts = self._qe_command.split()
            # Find and replace placeholders for input/output files
            for i, part in enumerate(command_parts):
                if input_filename in part:
                    command_parts[i] = part.replace(input_filename, str(input_filepath))
                if output_filename in part:
                    command_parts[i] = part.replace(output_filename, str(output_filepath))

            # The command should not include redirection, as we capture stdout
            command_parts = self._qe_command.split()
            # Find and replace placeholders for input/output files
            for i, part in enumerate(command_parts):
                if input_filename in part:
                    command_parts[i] = part.replace(input_filename, str(input_filepath))

            # Remove output file redirection if present
            command_parts = [p for p in command_parts if output_filename not in p and p != '>']

            process = subprocess.run(
                command_parts,
                check=True,
                cwd=str(temp_dir),
                capture_output=True,
                text=True,
            )

            # Write the captured output and then parse it
            output_filepath.write_text(process.stdout)
            result = parse_qe_output(process.stdout)

        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            # Handle cases where the QE run fails or output is not generated
            error_message = f"QE execution failed: {e}"
            if hasattr(e, "stderr"):
                error_message += f"\nStderr: {e.stderr}"
            output_str = (
                output_filepath.read_text() if output_filepath.exists() else ""
            )
            # Try parsing the output anyway, it might contain a QE error message
            result = parse_qe_output(output_str)
            if result.was_successful:  # If parser didn't find an error, add one
                result.was_successful = False
                result.error_message = error_message

        finally:
            # Clean up the temporary directory
            shutil.rmtree(temp_dir)

        # Write the final result to the database
        db_id = self._db.write(atoms, result)
        return db_id
