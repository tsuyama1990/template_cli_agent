import subprocess
import tempfile
import os
from typing import Dict
from ase.atoms import Atoms

from ..data.database import AseDB
from ..utils.dft_utils import generate_qe_input, parse_qe_output

class LabellingEngine:
    """
    Manages the execution of Quantum Espresso DFT calculations to label atomic structures.
    """
    def __init__(
        self,
        db: AseDB,
        qe_command: str,
        pseudo_dir: str,
        pseudopotentials: Dict[str, str],
        kpoints=(4, 4, 4),
        ecutwfc=60.0,
    ):
        """
        Initializes the LabellingEngine.

        Args:
            db: An instance of AseDB for database operations.
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x" or "mpirun -np 4 pw.x").
            pseudo_dir: Path to the directory containing pseudopotential files.
            pseudopotentials: A dictionary mapping element symbols to their pseudopotential filenames.
            kpoints: A tuple for the k-point grid.
            ecutwfc: The plane-wave cutoff energy in Rydberg.
        """
        self._db = db
        self._qe_command = qe_command
        self._pseudo_dir = pseudo_dir
        self._pseudopotentials = pseudopotentials
        self._kpoints = kpoints
        self._ecutwfc = ecutwfc

    def execute(self, atoms: Atoms) -> int:
        """
        Generates input, runs a QE calculation for the given Atoms object,
        parses the output, and saves the result to the database.

        Args:
            atoms: The ase.Atoms object to be calculated.

        Returns:
            The database ID of the new entry.
        """
        input_str = generate_qe_input(
            atoms=atoms,
            pseudo_dir=self._pseudo_dir,
            pseudopotentials=self._pseudopotentials,
            kpoints=self._kpoints,
            ecutwfc=self._ecutwfc,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            input_file_path = os.path.join(temp_dir, "qe_input.in")
            with open(input_file_path, "w") as f:
                f.write(input_str)

            command_list = self._qe_command.split() + ["-in", input_file_path]

            process = subprocess.run(
                command_list, capture_output=True, text=True, cwd=temp_dir
            )

            # Use stdout for parsing. If stdout is empty on error, fall back to stderr.
            output_to_parse = process.stdout
            if process.returncode != 0 and not output_to_parse:
                output_to_parse = process.stderr

            result = parse_qe_output(output_to_parse)

        # Write the result to the database
        db_id = self._db.write(atoms, result=result)
        return db_id
