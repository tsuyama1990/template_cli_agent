import io
import re
import subprocess

import ase.io
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.config.models import DFTParams


class QEOutputParsingError(Exception):
    """Custom exception for errors during QE output parsing."""
    pass


class LabellingEngine:
    """Encapsulates all interactions with Quantum Espresso."""

    def __init__(self, dft_config: DFTParams):
        """
        Initializes the LabellingEngine.

        Args:
            dft_config: The DFT-specific configuration.
        """
        self.dft_config = dft_config

    def run(self, atoms: Atoms) -> Atoms:
        """
        Runs a DFT calculation on a single Atoms object.

        Args:
            atoms: The structure to label.

        Returns:
            The same Atoms object with DFT results attached via a calculator.
        """
        input_str = self._generate_qe_input(atoms)
        success, stdout, stderr = self._execute_qe(input_str)

        if not success:
            # In a real implementation, we might try to recover here
            # or raise a more specific exception.
            raise RuntimeError(f"Quantum Espresso execution failed: {stderr}")

        results = self._parse_qe_output(stdout)

        # Attach results using a SinglePointCalculator
        calc = SinglePointCalculator(
            atoms,
            energy=results["energy"],
            forces=results["forces"],
            stress=results["stress"],
        )
        atoms.calc = calc
        return atoms

    def _generate_qe_input(self, atoms: Atoms) -> str:
        """
        Generates a QE input file string from an Atoms object.
        """
        # ASE's writer for QE requires specific parameters.
        input_params = {
            "control": {
                "calculation": "scf",
                "restart_mode": "from_scratch",
                "prefix": "pwscf",
                "outdir": "./out",
                "pseudo_dir": ".",  # Assumes pseudos are in the run dir
            },
            "system": {
                "ibrav": 0,
                "ecutwfc": self.dft_config.ecutwfc,
            },
            "electrons": {
                "mixing_beta": 0.7,
                "conv_thr": 1.0e-8,
            },
        }

        # We use io.StringIO to capture the output of ase.io.write as a string.
        with io.StringIO() as f:
            ase.io.write(
                f,
                atoms,
                format="espresso-in",
                pseudopotentials=self.dft_config.pseudopotentials,
                input_data=input_params,
            )
            return f.getvalue()

    def _execute_qe(self, input_str: str) -> tuple[bool, str, str]:
        """
        Executes the QE command with the given input string.
        """
        command = self.dft_config.command.split()
        try:
            process = subprocess.run(  # noqa: S603
                command,
                input=input_str,
                text=True,
                capture_output=True,
                check=False,  # We check the return code manually
                timeout=300,  # 5-minute timeout
            )
            success = process.returncode == 0
            return success, process.stdout, process.stderr
        except FileNotFoundError:
            msg = f"Command '{command[0]}' not found. Is QE installed and in PATH?"
            return False, "", msg
        except subprocess.TimeoutExpired:
            return False, "", "Quantum Espresso command timed out."

    def _parse_qe_output(self, output_text: str) -> dict:
        """
        Parses the QE output text to extract energy, forces, and stress.
        """
        # Conversion factor from kBar to eV/Angstrom^3
        KBAR_TO_EV_ANG3 = 1.0 / 160.21766208 * 0.1

        try:
            # --- Extract Energy ---
            energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output_text)
            if not energy_match:
                raise QEOutputParsingError("Final energy not found in QE output.")
            # Convert Ry to eV
            energy = float(energy_match.group(1)) * 13.605693122994

            # --- Extract Forces ---
            forces_block_match = re.search(
                r"Forces acting on atoms \(cartesian axes, Ry\/au\):\n\n([\s\S]+?)\n\n",
                output_text,
            )
            if not forces_block_match:
                raise QEOutputParsingError("Forces block not found in QE output.")

            forces_str = forces_block_match.group(1)
            forces = []
            for line in forces_str.strip().split("\n"):
                parts = line.split("=")
                force_components = [float(f) for f in parts[1].strip().split()]
                forces.append(force_components)
            # Convert Ry/au to eV/Angstrom
            forces = np.array(forces) * (13.605693122994 / 0.529177210903)

            # --- Extract Stress ---
            stress_pattern = (
                r"total   stress  \(Ry\/bohr\*\*3\)\s+\(kbar\)\n"
                r"([\s\S]+?)(?=\n\s*Pressure:)"
            )
            stress_block_match = re.search(stress_pattern, output_text)
            if not stress_block_match:
                raise QEOutputParsingError("Stress block not found in QE output.")

            stress_lines = stress_block_match.group(1).strip().split('\n')
            stress_matrix = np.array([list(map(float, line.split()[-3:])) for line in stress_lines])

            # Convert from kbar to eV/Angstrom^3
            stress_matrix *= KBAR_TO_EV_ANG3

            # Convert full 3x3 matrix to 6-element Voigt form for ASE
            stress_voigt = [
                stress_matrix[0, 0],
                stress_matrix[1, 1],
                stress_matrix[2, 2],
                stress_matrix[1, 2],
                stress_matrix[0, 2],
                stress_matrix[0, 1],
            ]

            return {
                "energy": energy,
                "forces": forces.tolist(),
                "stress": stress_voigt,
            }
        except (ValueError, IndexError) as e:
            raise QEOutputParsingError(f"Error parsing QE output: {e}") from e
