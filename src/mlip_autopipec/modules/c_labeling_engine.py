import os
import subprocess

from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.utils.dft_utils import create_qe_input_from_atoms, parse_qe_output


class LabelingEngine:
    """
    Manages the process of labeling atomic structures with DFT calculations.
    """

    def __init__(self, config: DFTCompute, db_wrapper: AseDBWrapper):
        """
        Initializes the LabelingEngine.

        Args:
            config: The DFTCompute configuration object.
            db_wrapper: An instance of the AseDBWrapper.
        """
        self.config = config
        self.db_wrapper = db_wrapper

    def execute(self):
        """
        Executes the labeling workflow for all unlabeled structures in the database.
        """
        print("Starting Labeling Engine...")
        rows_to_label = self.db_wrapper.get_rows_to_label()
        if not rows_to_label:
            print("No structures to label.")
            return

        print(f"Found {len(rows_to_label)} structures to label.")

        for i, row in enumerate(rows_to_label):
            atoms = row.toatoms()
            print(f"  Processing structure {i+1}/{len(rows_to_label)} (ID: {row.id})...")

            calc_dir = f"calc_{row.id}"
            os.makedirs(calc_dir, exist_ok=True)

            input_content = create_qe_input_from_atoms(
                atoms, self.config, self.config.pseudopotentials
            )
            input_file_path = os.path.join(calc_dir, 'qe.in')
            output_file_path = os.path.join(calc_dir, 'qe.out')

            with open(input_file_path, 'w') as f:
                f.write(input_content)

            command = self.config.command.split() + ['-in', input_file_path]

            try:
                with open(output_file_path, 'w') as out_f:
                    subprocess.run(  # noqa: S603
                        command,
                        stdout=out_f,
                        stderr=subprocess.PIPE,
                        check=True,
                        shell=False
                    )

                with open(output_file_path) as f:
                    output_content = f.read()

                dft_results = parse_qe_output(output_content)

                if dft_results:
                    print(
                        f"    ...DFT calculation successful. "
                        f"Energy: {dft_results['energy']:.4f} eV"
                    )
                    self.db_wrapper.update_row_with_dft_results(row.id, dft_results)
                else:
                    print("    ...DFT calculation failed: Could not parse output.")

            except FileNotFoundError:
                print(
                    f"    ...Error: Command '{self.config.command}' not found. "
                    f"Ensure Quantum Espresso is in your PATH."
                )
                break
            except subprocess.CalledProcessError as e:
                print(f"    ...DFT calculation failed with exit code {e.returncode}.")
                if e.stderr:
                    print(f"    ...Error output:\n{e.stderr.decode()}")
            finally:
                pass

        print("Labeling Engine finished.")
