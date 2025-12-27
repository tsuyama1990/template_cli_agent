"""The Labelling Engine (Module C) for automated DFT calculations."""
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils import dft_utils


class LabellingEngine:
    """
    Manages the execution of DFT calculations using Quantum Espresso.
    """

    def __init__(
        self,
        db: AseDB,
        qe_command: str,
        pseudo_dir: str,
        pseudos: Dict[str, str],
        kpts: Tuple[int, int, int],
        ecutwfc: float,
    ):
        """
        Initializes the Labelling Engine.

        Args:
            db: An instance of the AseDB wrapper.
            qe_command: The command to execute Quantum Espresso's pw.x.
                        (e.g., "pw.x" or "mpirun -np 4 pw.x").
            pseudo_dir: Path to the directory containing pseudopotential files.
            pseudos: A dictionary mapping atomic symbols to their pseudopotential filenames.
            kpts: A 3-tuple specifying the k-point grid for the SCF calculation.
            ecutwfc: The plane-wave cutoff energy in Rydberg.
        """
        self._db = db
        self._qe_command = qe_command
        self._pseudo_dir = pseudo_dir
        self._pseudos = pseudos
        self._kpts = kpts
        self._ecutwfc = ecutwfc

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation for a given atomic structure and saves the result.

        This method performs the following steps:
        1. Creates a temporary directory for the QE calculation.
        2. Generates the QE input file using the provided parameters.
        3. Executes the `pw.x` command using `subprocess`.
        4. Parses the stdout to extract energy, forces, and stress.
        5. Writes the result to the ASE database.

        Args:
            atoms: The ase.Atoms object representing the structure to be calculated.

        Returns:
            The unique database ID of the newly created record.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            work_dir = Path(temp_dir)
            input_file = work_dir / "qe.in"
            output_file = work_dir / "qe.out"

            # 1. Generate input file
            input_str = dft_utils.generate_qe_input(
                atoms=atoms,
                pseudo_dir=self._pseudo_dir,
                pseudos=self._pseudos,
                kpts=self._kpts,
                ecutwfc=self._ecutwfc,
                outdir=str(work_dir / "out"),
            )
            input_file.write_text(input_str)

            # 2. Execute pw.x
            command = self._qe_command.split() + ["-in", str(input_file)]
            process = subprocess.run(
                command,
                cwd=work_dir,
                capture_output=True,
                text=True,
            )
            output_file.write_text(process.stdout)

            # 3. Parse output
            if process.returncode != 0:
                result = dft_utils.DFTResult(
                    was_successful=False,
                    error_message=f"QE exited with error code {process.returncode}. Stderr: {process.stderr}"
                )
            else:
                result = dft_utils.parse_qe_output(process.stdout)

            # 4. Write to database
            db_id = self._db.write(atoms, result)
            return db_id
