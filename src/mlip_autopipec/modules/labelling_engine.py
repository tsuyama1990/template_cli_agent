"""The Labelling Engine for performing DFT calculations."""

import io

from ase import Atoms
from ase.io import write

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.utils import dft_utils
from mlip_autopipec.utils.runner import ProcessRunner


class LabellingEngine:
    """
    Manages the process of labelling atomic structures with DFT calculations.
    """

    def __init__(self, db: AseDB, runner: ProcessRunner, config: dict):
        """
        Initializes the LabellingEngine.

        Args:
            db: An AseDB instance for database interaction.
            runner: A ProcessRunner for executing external commands.
            config: A dictionary containing configuration for the engine.
        """
        self.db = db
        self.runner = runner
        self.config = config

    def run(self):
        """
        Runs the labelling process for all unlabelled structures in the database.
        """
        unlabelled_atoms = self.db.get_atoms_by_state("unlabelled")
        print(f"Found {len(unlabelled_atoms)} unlabelled structures to process.")

        for atoms in unlabelled_atoms:
            print(f"Processing structure with ID: {atoms.info['db_id']}")
            try:
                input_str = self._generate_qe_input(atoms)
                output_str = self._execute_dft(input_str)
                dft_result = dft_utils.parse_qe_output(output_str)
                self.db.update_with_dft_results(atoms.info["db_id"], dft_result)
                print(f"Successfully labelled structure ID: {atoms.info['db_id']}")
            except dft_utils.QEOutputError as e:
                print(f"Error processing structure ID {atoms.info['db_id']}: {e}")
            except Exception as e:
                print(f"An unexpected error occurred for structure ID {atoms.info['db_id']}: {e}")

    def _generate_qe_input(self, atoms: Atoms) -> str:
        """
        Generates a Quantum Espresso input file content as a string.

        NOTE: This is a simplified placeholder. A real implementation would have
        complex logic for selecting pseudopotentials, k-points, etc.
        """
        input_data = {
            "control": {
                "calculation": "scf",
                "prefix": "pwscf",
                "pseudo_dir": self.config.get("pseudo_dir", "."),
                "outdir": "./out",
            },
            "system": {
                "ibrav": 0,
                "nat": len(atoms),
                "ntyp": len(set(atoms.get_chemical_symbols())),
                "ecutwfc": self.config.get("ecutwfc", 50.0),
            },
            "electrons": {
                "mixing_beta": 0.7,
            },
        }

        string_io = io.StringIO()
        write(
            string_io,
            atoms,
            format="espresso-in",
            input_data=input_data,
            # Corrected argument name from pspots to pseudopotentials
            pseudopotentials={symbol: "pseudo.UPF" for symbol in set(atoms.get_chemical_symbols())},
            kpts=(1, 1, 1),
        )
        return string_io.getvalue()

    def _execute_dft(self, input_str: str) -> str:
        """
        Executes the DFT calculation.

        Args:
            input_str: The QE input file content.

        Returns:
            The standard output of the pw.x command.
        """
        command = self.config.get("dft_command", "pw.x").split()
        return self.runner.run(command, input_str)
