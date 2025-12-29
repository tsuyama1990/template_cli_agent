import subprocess
from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.utils import dft_utils


class LabellingEngine:
    def __init__(self, qe_command: str, db: AseDB):
        self._qe_command = qe_command.split() # Split command for shell=False
        self._db = db

    def execute(
        self,
        atoms: Atoms,
        parameters: dict[str, Any],
        pseudopotentials: dict[str, str],
        kpoints: tuple[int, int, int],
    ) -> int:
        """
        Generates input, runs a QE calculation, parses the output,
        and saves the result to the database.
        Returns the database ID of the new entry.
        """
        qe_input = dft_utils.generate_qe_input(
            atoms, parameters, pseudopotentials, kpoints
        )

        process_result = subprocess.run(
            self._qe_command,
            capture_output=True,
            text=True,
            input=qe_input,
        )

        dft_data = {}
        if process_result.returncode != 0:
            error_msg = (
                f"QE process failed with exit code {process_result.returncode}. "
                f"Stderr: {process_result.stderr}"
            )
            dft_data = {
                "was_successful": False,
                "error_message": error_msg,
            }
        else:
            parsed_output = dft_utils.parse_qe_output(process_result.stdout)
            if parsed_output["was_successful"]:
                dft_data = {
                    "total_energy_ev": parsed_output["total_energy_ev"],
                    "forces": parsed_output["forces"],
                    "stress": parsed_output["stress"],
                    "was_successful": True,
                }
            else:
                dft_data = {
                    "was_successful": False,
                    "error_message": parsed_output.get("error_message", "Unknown parsing error"),
                }

        # Create the Pydantic model first for validation
        dft_result = DFTResult.model_validate(dft_data)

        if dft_result.was_successful:
            voigt_stress = full_3x3_to_voigt_6_stress(np.array(dft_result.stress))
            calc = SinglePointCalculator(
                atoms,
                energy=dft_result.total_energy_ev,
                forces=np.array(dft_result.forces),
                stress=voigt_stress,
            )
            atoms.calc = calc

        db_id = self._db.write(atoms, dft_result)
        return db_id
