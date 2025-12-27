import subprocess
import shutil
from pathlib import Path

from ase.atoms import Atoms

from ..data.database import AseDB
from ..data.models import DFTResult
from ..utils.dft_utils import (
    generate_qe_input,
    parse_qe_output,
    SSSP_PSEUDOS,
    SSSP_CUTOFFS
)

class LabellingEngine:
    """
    The core engine for performing automated DFT calculations (labelling).
    This class orchestrates the process of taking an atomic structure,
    running a Quantum Espresso calculation, and storing the results.
    """

    def __init__(self, qe_command: str, db: AseDB):
        """
        Initialises the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso's pw.x.
                        This should include any MPI prefixes, e.g., "mpirun -np 4 pw.x".
            db: An instance of the AseDB wrapper class for database access.
        """
        if not shutil.which(qe_command.split()[0]):
            raise FileNotFoundError(
                f"The QE executable '{qe_command.split()[0]}' was not found in the system's PATH. "
                "Please ensure Quantum Espresso is installed and accessible."
            )
        self._qe_command = qe_command
        self._db = db

    def execute(self, atoms: Atoms, work_dir: Path = Path("./qe_temp")) -> int:
        """
        Executes the full labelling workflow for a given atomic structure.

        This method performs the following steps:
        1. Creates a temporary working directory for the QE calculation.
        2. Determines the appropriate pseudopotentials and cutoffs.
        3. Generates the QE input file using the utility function.
        4. Runs the QE calculation as a subprocess.
        5. Parses the output to extract energy, forces, and stress.
        6. Populates a DFTResult data model with the results.
        7. Writes the structure and results to the database.
        8. Cleans up the temporary directory.

        Args:
            atoms: The ase.Atoms object to be calculated.
            work_dir: The directory to perform the calculation in.

        Returns:
            The integer ID of the new row in the database.
        """
        work_dir.mkdir(exist_ok=True)
        input_file = work_dir / "qe_input.in"
        output_file = work_dir / "qe_output.out"

        try:
            # For simplicity, we use predefined SSSP pseudos and cutoffs.
            # A more advanced implementation would fetch these from a config.
            symbols = set(atoms.get_chemical_symbols())
            pseudos = {s: SSSP_PSEUDOS[s] for s in symbols}

            # Use the maximum cutoff required by any element present
            ecutwfc = max(SSSP_CUTOFFS[s][0] for s in symbols)
            ecutrho = max(SSSP_CUTOFFS[s][1] for s in symbols)

            # A simple k-point grid selection heuristic
            k_points = (4, 4, 4) if atoms.pbc.all() else (1, 1, 1)

            input_str = generate_qe_input(atoms, pseudos, ecutwfc, ecutrho, k_points)
            input_file.write_text(input_str)

            # Construct the full command for the subprocess
            command = self._qe_command.split() + ["-in", str(input_file)]

            # Execute the QE calculation
            process_result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                cwd=work_dir
            )

            output_file.write_text(process_result.stdout)

            # Parse the results from stdout
            energy, forces, stress, success, parsed_error = parse_qe_output(process_result.stdout)

            final_error = parsed_error
            # If the subprocess failed but the parser found no specific error,
            # create a generic one.
            if process_result.returncode != 0 and not final_error:
                final_error = f"QE process exited with non-zero status {process_result.returncode}."

            # A non-zero return code or a parsed error message means the run failed.
            if process_result.returncode != 0 or final_error:
                success = False

            dft_result = DFTResult(
                total_energy_ev=energy,
                forces=forces,
                stress=stress,
                was_successful=success,
                error_message=final_error
            )

        except Exception as e:
            # Catch any other unexpected errors during the process
            dft_result = DFTResult(
                total_energy_ev=None,
                forces=None,
                stress=None,
                was_successful=False,
                error_message=str(e)
            )
        finally:
            # Ensure cleanup happens even if there's an error
            # For debugging, one might want to comment this line out
            shutil.rmtree(work_dir)

        # Write the final result to the database
        db_id = self._db.write(atoms, result=dft_result)
        return db_id
