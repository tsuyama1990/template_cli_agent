import os
import subprocess
import tempfile
from ase.calculators.espresso import Espresso
from ase.io import read
import numpy as np

from ..database import AseDBWrapper
from ..config import DFTInputConfig, DFTResult

class LabelingEngine:
    """Handles the execution of DFT calculations for labeling."""

    def __init__(self, db_wrapper: AseDBWrapper, dft_config: DFTInputConfig, qe_command: str):
        """
        Initializes the LabelingEngine.

        Args:
            db_wrapper: An instance of the AseDBWrapper.
            dft_config: Configuration for the DFT calculation.
            qe_command: The command to execute Quantum Espresso's pw.x.
        """
        self.db_wrapper = db_wrapper
        self.dft_config = dft_config
        self.qe_command = qe_command

    def label_structure(self, uid: int):
        """
        Labels a single atomic structure in the database.

        Args:
            uid: The ID of the structure to label.
        """
        print(f"Starting labeling for structure ID: {uid}")
        atoms = self.db_wrapper.get_atoms_by_id(uid)
        if not atoms:
            print(f"Error: Structure with ID {uid} not found.")
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            pseudopotentials = self.dft_config.pseudopotentials

            # ASE's Espresso calculator requires the command to be set via an env var
            # or passed to the constructor. We'll use the constructor.
            # It also needs the `label` parameter for file naming.
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

            try:
                # This call triggers the calculation
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces()
                # Get the full 3x3 stress tensor
                stress = atoms.get_stress(voigt=False)

                dft_result = DFTResult(
                    energy=energy,
                    forces=forces,
                    stress=stress
                )

                self.db_wrapper.update_labels(uid, dft_result)
                print(f"Successfully labeled structure ID: {uid}")
                print(f"  - Energy: {energy:.4f} eV")

            except subprocess.CalledProcessError as e:
                print(f"Error during Quantum Espresso execution for ID {uid}:")
                print(e.stdout)
                print(e.stderr)
                # Optionally, update the state to 'labeling_failed'
                # self.db_wrapper.update_state(uid, "labeling_failed")
            except Exception as e:
                print(f"An unexpected error occurred during labeling of ID {uid}: {e}")
