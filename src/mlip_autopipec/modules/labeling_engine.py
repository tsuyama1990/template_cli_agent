import os
import subprocess
import tempfile
import logging
from ase.calculators.espresso import Espresso
from ase.io import read
import numpy as np

from ..database import AseDBWrapper
from ..config import DFTInputConfig, DFTResult

# Set up a logger for this module
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LabelingEngine:
    """
    Handles the execution of DFT calculations for labeling atomic structures.

    This engine is responsible for taking a structure from the database,
    running a Quantum Espresso calculation using ASE's built-in calculator,
    parsing the results, and storing them back into the database.
    """

    def __init__(self, db_wrapper: AseDBWrapper, dft_config: DFTInputConfig, qe_command: str):
        """
        Initializes the LabelingEngine.

        Args:
            db_wrapper: An instance of the AseDBWrapper for database interactions.
            dft_config: A Pydantic model containing the configuration for the
                        DFT calculation (e.g., k-points, cutoffs).
            qe_command: The shell command used to execute the Quantum Espresso
                        `pw.x` binary (e.g., 'pw.x -in input.pwi > output.pwo').
        """
        self.db_wrapper = db_wrapper
        self.dft_config = dft_config
        self.qe_command = qe_command

    def label_structure(self, uid: int):
        """
        Labels a single atomic structure identified by its database ID.

        This method performs the following steps:
        1. Retrieves the atomic structure from the database.
        2. Creates a temporary directory for the DFT calculation.
        3. Configures and attaches an ASE `Espresso` calculator.
        4. Executes the DFT calculation via a subprocess call.
        5. Parses the resulting energy, forces, and stress.
        6. Populates a `DFTResult` object.
        7. Updates the structure's entry in the database with the results.

        If the calculation fails, it logs the error but does not currently
        update the database state to 'labeling_failed'.

        Args:
            uid: The unique ID of the structure in the database to be labeled.
        """
        logger.info(f"Starting labeling for structure ID: {uid}")
        atoms = self.db_wrapper.get_atoms_by_id(uid)
        if not atoms:
            logger.error(f"Structure with ID {uid} not found.")
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            pseudopotentials = self.dft_config.pseudopotentials

            calc = Espresso(
                command=self.qe_command,
                label="qe_calc",
                pseudopotentials=pseudopotentials,
                kpts=self.dft_config.kpoints,
                ecutwfc=self.dft_config.ecutwfc,
                input_data=self.dft_config.control,
                directory=tmpdir
            )

            atoms.calc = calc

            logger.info(f"Executing Quantum Espresso command for ID {uid}: {self.qe_command}")
            try:
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces()
                stress = atoms.get_stress(voigt=False)

                dft_result = DFTResult(
                    energy=energy,
                    forces=forces,
                    stress=stress
                )

                self.db_wrapper.update_labels(uid, dft_result)
                logger.info(f"Successfully labeled structure ID: {uid} with Energy: {energy:.4f} eV")

            except subprocess.CalledProcessError as e:
                logger.error(f"Error during Quantum Espresso execution for ID {uid}:")
                logger.error(f"STDOUT: {e.stdout}")
                logger.error(f"STDERR: {e.stderr}")
            except Exception as e:
                logger.error(f"An unexpected error occurred during labeling of ID {uid}: {e}", exc_info=True)
