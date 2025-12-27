import subprocess
import tempfile
from pathlib import Path
from typing import Any

from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output


class LabellingEngine:
    """The engine for performing DFT labelling."""

    def __init__(
        self,
        qe_command: str,
        db: AseDB,
        parameters: dict[str, Any],
        pseudopotentials: dict[str, str],
    ):
        """
        Initializes the LabellingEngine.

        Args:
            qe_command: The command to execute Quantum Espresso (e.g., "pw.x").
            db: An instance of AseDB.
            parameters: A dictionary of QE parameters.
            pseudopotentials: A mapping from atomic symbol to pseudopotential filename.
        """
        self._qe_command = qe_command
        self._db = db
        self._parameters = parameters
        self._pseudopotentials = pseudopotentials

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation for the given Atoms object and saves the result to the database.

        Args:
            atoms: The ASE Atoms object to label.

        Returns:
            The database ID of the new entry.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_file = tmp_path / "qe.in"
            output_file = tmp_path / "qe.out"

            input_content = generate_qe_input(atoms, self._parameters, self._pseudopotentials)
            input_file.write_text(input_content)

            command = f"{self._qe_command} -in {input_file} > {output_file}"

            try:
                subprocess.run(command, shell=True, check=True, cwd=tmp_path)
                output_content = output_file.read_text()
                was_successful, error_message, results = parse_qe_output(output_content)

                if was_successful:
                    dft_result = DFTResult(
                        total_energy_ev=results["total_energy_ev"],
                        forces=results["forces"],
                        stress=results["stress"],
                        was_successful=True,
                    )
                else:
                    dft_result = DFTResult(
                        total_energy_ev=0.0,
                        forces=[[0.0, 0.0, 0.0]] * len(atoms),
                        stress=[[0.0, 0.0, 0.0]] * 3,
                        was_successful=False,
                        error_message=error_message,
                    )

            except subprocess.CalledProcessError as e:
                dft_result = DFTResult(
                    total_energy_ev=0.0,
                    forces=[[0.0, 0.0, 0.0]] * len(atoms),
                    stress=[[0.0, 0.0, 0.0]] * 3,
                    was_successful=False,
                    error_message=f"Quantum Espresso execution failed: {e}",
                )

        db_id = self._db.write(atoms, dft_result)
        return db_id
