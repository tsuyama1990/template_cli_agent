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

            calc_dir = f"calc_{row_id}"
            os.makedirs(calc_dir, exist_ok=True)

            input_content = create_qe_input_from_atoms(
                atoms, self.config, self.config.pseudopotentials
            )
            input_file_path = os.path.join(calc_dir, QE_INPUT_FILENAME)
            output_file_path = os.path.join(calc_dir, QE_OUTPUT_FILENAME)

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
                        shell=False,
                    )

                with open(output_file_path) as f:
                    output_content = f.read()

                dft_results = parse_qe_output(output_content)

                if dft_results:
                    logger.info(
                        f"    ...DFT calculation successful. "
                        f"Energy: {dft_results.energy:.4f} eV"
                    )
                    results.append((row_id, dft_results))
                else:
                    logger.warning(
                        f"    ...DFT calculation for ID {row_id} failed: "
                        "Could not parse output."
                    )

            except FileNotFoundError as e:
                logger.error(
                    f"    ...Error: Command '{self.config.command}' not found. "
                    "Ensure Quantum Espresso is in your PATH."
                )
                raise FileNotFoundError(
                    f"Command '{self.config.command}' not found."
                ) from e
            except subprocess.CalledProcessError as e:
                logger.error(
                    f"    ...DFT calculation for ID {row_id} failed with exit code "
                    f"{e.returncode}."
                )
                if e.stderr:
                    logger.error(f"    ...Error output:\n{e.stderr.decode()}")
                raise subprocess.CalledProcessError(
                    e.returncode, e.cmd, output=e.output, stderr=e.stderr
                ) from e

        logger.info("Labeling Engine finished.")
        return results
