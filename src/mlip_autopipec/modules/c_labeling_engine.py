import subprocess
import tempfile
import os
import ase.io
from ase.calculators.espresso import Espresso

from mlip_autopipec.configs.models import DFTComputeConfig
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.utils.dft_utils import parse_qe_output


class LabelingEngine:
    """
    Manages the process of labeling atomic structures with DFT calculations.
    """

    def __init__(self, config: DFTComputeConfig, db_wrapper: AseDBWrapper):
        self.config = config
        self.db = db_wrapper

    def run(self):
        """
        Fetches unlabeled structures, runs DFT calculations, and updates the database.
        """
        print("Starting Labeling Engine...")
        unlabeled_rows = self.db.get_unlabeled_rows()
        print(f"Found {len(unlabeled_rows)} unlabeled structures to process.")

        for row in unlabeled_rows:
            atoms = row.toatoms()
            print(f"Processing structure with ID: {row.id}")

            with tempfile.TemporaryDirectory() as tmpdir:
                input_path = os.path.join(tmpdir, "qe_input.in")
                output_path = os.path.join(tmpdir, "qe_output.out")

                pseudopotentials = {
                    symbol: f"{symbol}.{self.config.pseudopotentials}.UPF"
                    for symbol in set(atoms.get_chemical_symbols())
                }

                input_data = {
                    "control": {
                        "pseudo_dir": os.environ.get("PSEUDO_DIR", "."),
                        "calculation": "scf",
                        "tstress": True,
                        "tprnfor": True,
                    },
                    "system": {
                        "ecutwfc": self.config.ecutwfc,
                        "ecutrho": self.config.ecutrho,
                        "smearing": self.config.smearing,
                        "degauss": self.config.degauss,
                    },
                    "electrons": {
                        "mixing_beta": 0.7
                    }
                }

                # Create and configure the Espresso calculator
                calc = Espresso(
                    pseudopotentials=pseudopotentials,
                    kpts=None,  # Let ASE generate k-points from density
                    kdensity=self.config.kpoints_density,
                    input_data=input_data,
                    label="espresso" # ASE default prefix
                )

                # Attach the calculator to the atoms object
                atoms.calc = calc

                # Use ase.io.write to generate the input file from the attached calculator
                ase.io.write(input_path, atoms, format="espresso-in")

                # Execute Quantum Espresso
                command = self.config.command.split() + ["-in", input_path]
                with open(output_path, "w") as f_out:
                    result = subprocess.run(
                        command,
                        stdout=f_out,
                        stderr=subprocess.PIPE,
                        text=True,
                        cwd=tmpdir,
                    )

                if result.returncode != 0:
                    print(f"Error running QE for structure {row.id}:")
                    print(result.stderr)
                    continue

                # Read the output file for parsing
                with open(output_path, "r") as f_out:
                    output_content = f_out.read()

                # Parse the output
                try:
                    parsed_data = parse_qe_output(output_content)
                except ValueError as e:
                    print(f"Error parsing QE output for structure {row.id}: {e}")
                    continue

                # Update the database
                self.db.update_row(
                    row_id=row.id,
                    data=parsed_data,
                    key_value_pairs={"state": "labeled"},
                )
                print(f"Successfully labeled structure {row.id}.")

        print("Labeling Complete.")
