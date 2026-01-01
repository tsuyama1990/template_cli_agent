import logging
import os
import subprocess

from ase import Atoms

from mlip_autopipec.data.models import DFTCompute, DFTResults
from mlip_autopipec.utils.dft_utils import create_qe_input_from_atoms, parse_qe_output

QE_INPUT_FILENAME = "qe.in"
QE_OUTPUT_FILENAME = "qe.out"

logger = logging.getLogger(__name__)


class LabelingEngine:
    """
    Manages the process of labeling atomic structures with DFT calculations.
    """

    def __init__(self, config: DFTCompute):
        """
        Initializes the LabelingEngine.

        Args:
            config: The DFTCompute configuration object.
        """
        self.config = config

    def execute(
        self, structures: list[tuple[int, Atoms]]
    ) -> list[tuple[int, DFTResults]]:
        """
        Executes the labeling workflow for the given structures.

        Args:
            structures: A list of tuples, each containing a row ID and an ASE Atoms object.

        Returns:
            A list of tuples, each containing a row ID and a DFTResults object.
        """
        logger.info("Starting Labeling Engine...")
        if not structures:
            logger.info("No structures to label.")
            return []

        logger.info(f"Received {len(structures)} structures to label.")
        results = []
        for i, (row_id, atoms) in enumerate(structures):
            logger.info(
                f"  Processing structure {i+1}/{len(structures)} (ID: {row_id})..."
            )
            try:
                dft_results = self._run_single_calculation(row_id, atoms)
                if dft_results:
                    results.append((row_id, dft_results))
            except (OSError, subprocess.SubprocessError) as e:
                logger.error(f"    ...Critical error processing structure ID {row_id}. "
                             f"Skipping this structure. Error: {e}", exc_info=True)
                # Continue to the next structure
                continue

        logger.info("Labeling Engine finished.")
        return results

    def _run_single_calculation(self, row_id: int, atoms: Atoms) -> DFTResults | None:
        """Runs a single DFT calculation for a given atoms object."""
        calc_dir = f"calc_{row_id}"
        try:
            os.makedirs(calc_dir, exist_ok=True)
        except OSError as e:
            logger.error(f"    ...Could not create calculation directory '{calc_dir}'.")
            raise OSError(f"Failed to create directory '{calc_dir}'") from e

        input_file_path = os.path.join(calc_dir, QE_INPUT_FILENAME)
        output_file_path = os.path.join(calc_dir, QE_OUTPUT_FILENAME)

        try:
            input_content = create_qe_input_from_atoms(
                atoms, self.config, self.config.pseudopotentials
            )
            with open(input_file_path, 'w') as f:
                f.write(input_content)
        except OSError as e:
            logger.error(f"    ...Could not write to input file '{input_file_path}'.")
            raise OSError(f"Failed to write file '{input_file_path}'") from e

        command = self.config.command.split() + ['-in', input_file_path]

        try:
            with open(output_file_path, 'w') as out_f:
                subprocess.run(  # noqa: S603
                    command,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    check=True,
                    shell=False,
                )
        except FileNotFoundError as e:
            logger.error(
                f"    ...Error: Command '{self.config.command}' not found. "
                "Ensure Quantum Espresso is in your PATH."
            )
            raise FileNotFoundError(f"Command '{self.config.command}' not found.") from e
        except subprocess.CalledProcessError as e:
            logger.error(
                f"    ...DFT calculation for ID {row_id} failed with exit code "
                f"{e.returncode}."
            )
            if e.stderr:
                logger.error(f"    ...Error output:\n{e.stderr.decode()}")
            raise e # Re-raise to be caught by the main loop

        try:
            with open(output_file_path) as f:
                output_content = f.read()
        except OSError as e:
            logger.error(f"    ...Could not read output file '{output_file_path}'.")
            raise OSError(f"Failed to read file '{output_file_path}'") from e

        dft_results = parse_qe_output(output_content)
        if dft_results:
            logger.info(
                f"    ...DFT calculation successful. "
                f"Energy: {dft_results.energy:.4f} eV"
            )
            return dft_results

        logger.warning(
            f"    ...DFT calculation for ID {row_id} finished, "
            "but output could not be parsed."
        )
        return None
