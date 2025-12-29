# ruff: noqa: D101, D102, D103, D104, D105, D107
import subprocess
import tempfile
from pathlib import Path

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output


class LabellingEngine:
    """Handles the execution of DFT calculations using Quantum Espresso."""

    def __init__(
        self,
        qe_command: str,
        db: AseDB,
        pseudo_potentials: dict[str, str],
        k_points: tuple[int, int, int] = (1, 1, 1),
        ecutwfc: float = 60.0,
    ):
        """Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x").
            db: An instance of the AseDB wrapper.
            pseudo_potentials: A dictionary mapping element symbols to pseudopotential filenames.
            k_points: A tuple of 3 integers for the k-point grid.
            ecutwfc: The plane-wave cutoff energy in Ry.
        """
        self._qe_command = (
            qe_command.split() if isinstance(qe_command, str) else qe_command
        )
        self._db = db
        self._pseudo_potentials = pseudo_potentials
        self._k_points = k_points
        self._ecutwfc = ecutwfc

    def execute(self, atoms: Atoms) -> int:
        """
        Generates input, runs QE, parses output, and saves the result to the DB.
        Returns the database ID of the new entry.
        """
        input_str = generate_qe_input(
            atoms, self._pseudo_potentials, self._k_points, self._ecutwfc
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_file = temp_path / "qe.in"
            output_file = temp_path / "qe.out"
            (temp_path / "out").mkdir()  # QE 'outdir'

            with open(input_file, "w") as f:
                f.write(input_str)

            # Construct the command for subprocess
            command = self._qe_command + ["-in", str(input_file)]

            # Run the QE calculation
            process_result = subprocess.run(  # nosec
                command,
                capture_output=True,
                text=True,
                cwd=temp_path,
            )
            with open(output_file, "w") as f:
                f.write(process_result.stdout)

            output_str = process_result.stdout
            result = parse_qe_output(output_str)

        db_id = self._db.write(atoms, result)
        return db_id
