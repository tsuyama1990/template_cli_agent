import subprocess
import tempfile
import os
from typing import Dict, Any

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

class LabellingEngine:
    def __init__(self, qe_command: str, parameters: Dict[str, Any], db: AseDB):
        self.qe_command = qe_command
        self.parameters = parameters
        self.db = db

    def execute(self, atoms: Atoms) -> int:
        """
        Runs a DFT calculation, parses the output, and stores the result in the database.
        Returns the database ID of the new entry.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, 'qe.in')
            input_str = generate_qe_input(atoms, self.parameters)
            with open(input_path, 'w') as f:
                f.write(input_str)

            command = f"{self.qe_command} -in {input_path}"
            was_successful = False
            error_message = None

            try:
                process = subprocess.run(
                    command, shell=True, capture_output=True, text=True, check=True, cwd=temp_dir
                )
                output_str = process.stdout
                parsed_data, error = parse_qe_output(output_str)
                if error:
                    error_message = error
                else:
                    was_successful = True
            except subprocess.CalledProcessError as e:
                output_str = (e.stdout or "") + (e.stderr or "")
                error_message = f"QE execution failed. Output:\n{output_str}"

            metadata = {'was_successful': was_successful, 'error_message': error_message}

            if was_successful and parsed_data:
                calc = SinglePointCalculator(
                    atoms,
                    energy=parsed_data['total_energy_ev'],
                    forces=parsed_data['forces'],
                    stress=parsed_data['stress']
                )
                atoms.calc = calc

            db_id = self.db.write(atoms, metadata=metadata)

        return db_id
