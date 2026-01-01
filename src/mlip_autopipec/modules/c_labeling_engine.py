# src/mlip_autopipec/modules/c_labeling_engine.py

import io
import subprocess
from pathlib import Path

from ase.io import write

from mlip_autopipec.configs.models import DFTComputeConfig
from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.utils.dft_utils import parse_qe_output


class LabelingEngine:
    """An engine for labeling atomic structures with DFT calculations."""

    def __init__(
        self,
        config: DFTComputeConfig,
        db_wrapper: AseDBWrapper,
        calculation_dir: Path = Path("calculations"),
    ):
        """
        Initializes the LabelingEngine.

        Args:
            config: The DFT compute configuration.
            db_wrapper: The database wrapper instance.
            calculation_dir: The directory to store calculation files.
        """
        self.config = config
        self.db_wrapper = db_wrapper
        self.calculation_dir = calculation_dir
        self.calculation_dir.mkdir(exist_ok=True)

    def run(self):
        """
        Runs the labeling process for all unlabeled structures in the database.
        """
        unlabeled_rows = self.db_wrapper.get_unlabeled_rows()
        print(f"Found {len(unlabeled_rows)} unlabeled structures to process.")

        for row in unlabeled_rows:
            atoms = row.toatoms()
            row_id = row.id
            print(f"Processing structure with ID: {row_id}")

            # Create a dedicated directory for this calculation
            single_calc_dir = self.calculation_dir / f"id_{row_id}"
            single_calc_dir.mkdir(exist_ok=True)
            input_path = single_calc_dir / "pw.in"
            output_path = single_calc_dir / "pw.out"

            # Generate QE input file
            pseudopotentials = {atom.symbol: "" for atom in atoms}
            kpts = (
                int(atoms.cell.lengths()[0] * self.config.kpoints_density),
                int(atoms.cell.lengths()[1] * self.config.kpoints_density),
                int(atoms.cell.lengths()[2] * self.config.kpoints_density),
            )

            input_data = {
                "control": {"calculation": "vc-relax"},
                "system": {
                    "ecutwfc": self.config.ecutwfc,
                    "ecutrho": self.config.ecutrho,
                    "occupations": "smearing",
                    "smearing": self.config.smearing,
                    "degauss": self.config.degauss,
                },
            }

            with io.StringIO() as buffer:
                write(
                    buffer,
                    atoms,
                    format="espresso-in",
                    pseudopotentials=pseudopotentials,
                    kpts=kpts,
                    input_data=input_data,
                )
                input_file_content = buffer.getvalue()

            with open(input_path, "w") as f:
                f.write(input_file_content)

            # Run Quantum Espresso
            command = self.config.command.split() + ["-in", str(input_path)]
            try:
                result = subprocess.run(  # noqa: S603
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                    cwd=single_calc_dir,
                )
                with open(output_path, "w") as f:
                    f.write(result.stdout)

                # Parse output and update database
                parsed_results = parse_qe_output(result.stdout)
                if parsed_results:
                    # ASE's DB requires results in a 'data' dictionary
                    # And metadata in key_value_pairs
                    dft_data = {
                        "energy": parsed_results.get("energy"),
                        "forces": parsed_results.get("forces"),
                        "stress": parsed_results.get("stress"),
                    }
                    metadata = {"state": "labeled"}
                    self.db_wrapper.update_row(row_id, dft_data, metadata)
                    print(f"Successfully labeled structure {row_id}.")
                else:
                    print(f"Error: Failed to parse QE output for structure {row_id}.")
                    metadata = {"state": "error"}
                    self.db_wrapper.update_row(row_id, {}, metadata)

            except subprocess.CalledProcessError as e:
                print(f"Error running QE for structure {row_id}: {e}")
                with open(output_path, "w") as f:
                    f.write(e.stdout)
                    f.write("\n--- STDERR ---\n")
                    f.write(e.stderr)
                metadata = {"state": "error"}
                self.db_wrapper.update_row(row_id, {}, metadata)
