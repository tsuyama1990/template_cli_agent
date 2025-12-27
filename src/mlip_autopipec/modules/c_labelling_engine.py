import subprocess
import tempfile
import os
from ase.atoms import Atoms
from typing import Dict, Any

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

class LabellingEngine:
    def __init__(self, qe_command: str, db: AseDB, default_params: Dict[str, Any], pseudos: Dict[str, str]):
        self._qe_command = qe_command
        self._db = db
        self._default_params = default_params
        self._pseudos = pseudos

    def execute(self, atoms: Atoms) -> int:
        """
        Generates input, runs QE, parses output, and saves the result to the database.
        Returns the database ID of the new entry.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, 'qe.in')
            output_file = os.path.join(temp_dir, 'qe.out')

            # Add temp_dir to control params for QE
            params = self._default_params.copy()
            params['control'] = params.get('control', {})
            params['control']['outdir'] = temp_dir
            params['control']['pseudo_dir'] = '.' # Assuming pseudos are in the working dir or path

            input_str = generate_qe_input(atoms, params, self._pseudos)
            with open(input_file, 'w') as f:
                f.write(input_str)

            # Construct and run the QE command
            command = f"{self._qe_command} -in {input_file}"
            try:
                with open(output_file, 'w') as f_out:
                    process = subprocess.run(
                        command,
                        shell=True,
                        stdout=f_out,
                        stderr=subprocess.PIPE,
                        check=True,
                        text=True,
                    )
                with open(output_file, 'r') as f_out:
                    output_str = f_out.read()
                dft_result = parse_qe_output(output_str)

            except subprocess.CalledProcessError as e:
                with open(output_file, 'r') as f_out:
                    output_str = f_out.read()
                dft_result = parse_qe_output(output_str) # Still try to parse for QE's error message
                if not dft_result.error_message: # If parser didn't find a specific QE error
                     dft_result = DFTResult(
                        total_energy_ev=0.0,
                        forces=[],
                        stress=[],
                        was_successful=False,
                        error_message=f"QE execution failed with exit code {e.returncode}. Stderr: {e.stderr}"
                     )

        # Write the result to the database regardless of success
        db_id = self._db.write(atoms, dft_result)
        return db_id
