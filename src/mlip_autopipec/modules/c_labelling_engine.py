from ase.atoms import Atoms
import subprocess
from typing import Dict, Tuple

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output
from mlip_autopipec.data.models import DFTResult

class LabellingEngine:
    """
    Manages the execution of DFT calculations using Quantum Espresso.
    """
    def __init__(self, qe_command: str, db: AseDB, pseudos: Dict[str, str], kpts: Tuple[int, int, int], ecutwfc: float):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "mpirun -np 4 pw.x").
            db: An instance of the AseDB class.
            pseudos: A dictionary mapping element symbols to pseudopotential filenames.
            kpts: A tuple for the k-point grid (e.g., (4, 4, 4)).
            ecutwfc: The plane-wave cutoff energy in Ry.
        """
        self._qe_command = qe_command
        self._db = db
        self._pseudos = pseudos
        self._kpts = kpts
        self._ecutwfc = ecutwfc

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation for the given atoms object, parses the output,
        and saves the result to the database.

        Args:
            atoms: The ASE Atoms object to be calculated.

        Returns:
            The database ID of the new entry.
        """
        # 1. Generate Quantum Espresso input file content
        input_data = generate_qe_input(
            atoms=atoms,
            pseudos=self._pseudos,
            kpts=self._kpts,
            ecutwfc=self._ecutwfc,
        )

        # 2. Execute Quantum Espresso
        # In a real scenario, this would involve managing files and directories.
        # For simplicity, we pass input via stdin and capture stdout.
        try:
            process = subprocess.run(
                self._qe_command.split(),
                input=input_data,
                capture_output=True,
                text=True,
                check=True,
                timeout=600, # 10 minute timeout
            )
            stdout = process.stdout

        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
            # Handle cases where pw.x fails, is not found, or times out
            result = DFTResult(
                total_energy_ev=0.0, forces=[], stress=[], was_successful=False,
                error_message=f"Subprocess execution failed: {e}"
            )
            return self._db.write(atoms, result)

        # 3. Parse the output
        result = parse_qe_output(stdout)

        # 4. Save to the database
        db_id = self._db.write(atoms, result)
        return db_id
