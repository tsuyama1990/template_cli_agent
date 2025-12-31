import subprocess
from pathlib import Path

from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.utils.dft_utils import create_qe_input_from_atoms, parse_qe_output


class LabelingEngine:
    """Engine for labeling atomic structures with DFT calculations."""

    def __init__(self, db_wrapper: AseDBWrapper, config: DFTCompute):
        self.db_wrapper = db_wrapper
        self.config = config

    def execute(self):
        """Executes the labeling workflow for all unlabeled structures."""
        rows_to_label = self.db_wrapper.get_rows_to_label()
        for row in rows_to_label:
            self._run_dft_calculation(row)

    def _run_dft_calculation(self, row):
        """Runs a single DFT calculation for a given database row."""
        atoms = row.toatoms()
        calc_dir = Path(f"calc_{row.id}")
        calc_dir.mkdir(exist_ok=True)

        input_file = calc_dir / "pw.in"
        output_file = calc_dir / "pw.out"

        input_str = create_qe_input_from_atoms(atoms, self.config)
        input_file.write_text(input_str)

        command = self.config.command.split() + ["-in", str(input_file)]

        with open(output_file, "w") as f_out:
            subprocess.run(
                command,
                stdout=f_out,
                stderr=subprocess.PIPE,
                check=True,
                cwd=calc_dir,
            )

        output_str = output_file.read_text()
        dft_results = parse_qe_output(output_str)

        if "energy" in dft_results and "forces" in dft_results:
            calc = SinglePointCalculator(
                atoms,
                energy=dft_results.get("energy"),
                forces=dft_results.get("forces"),
                stress=dft_results.get("stress"),
            )
            atoms.calc = calc
            self.db_wrapper.update_row_with_dft_results(row.id, atoms)
