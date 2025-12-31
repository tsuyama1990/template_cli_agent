"""Utility functions for DFT calculations, particularly for Quantum Espresso."""

import re
import numpy as np
from mlip_autopipec.domain.models import DFTResult


class QEOutputError(Exception):
    """Custom exception for errors during QE output parsing."""
    pass


def parse_qe_output(output_str: str) -> DFTResult:
    """
    Parses the text output of a Quantum Espresso pw.x calculation.

    Args:
        output_str: The stdout content from the pw.x run.

    Returns:
        A DFTResult object containing the extracted energy, forces, and stress.

    Raises:
        QEOutputError: If parsing fails or the calculation did not converge.
    """
    if "convergence NOT achieved" in output_str:
        raise QEOutputError("SCF calculation did not converge.")

    energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output_str)
    if not energy_match:
        raise QEOutputError("Could not find total energy in QE output.")
    energy = float(energy_match.group(1))

    forces_pattern = (
        r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n\s*\n"
        r"((?:\s*atom\s+\d+\s+force\s+=\s+[-.\d\s]+\n)+)"
    )
    forces_match = re.search(forces_pattern, output_str)
    if not forces_match:
        raise QEOutputError("Could not find forces in QE output.")

    forces_lines = forces_match.group(1).strip().split('\n')
    forces = []
    for line in forces_lines:
        parts = line.split()
        forces.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])

    stress_pattern = (
        r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*[-.\d]+\n"
        r"((?:\s*[-.\d\s]+\n){3})"
    )
    stress_match = re.search(stress_pattern, output_str)
    if not stress_match:
        raise QEOutputError("Could not find stress tensor in QE output.")

    stress_lines = stress_match.group(1).strip().split('\n')
    stress = []
    for line in stress_lines:
        stress.append([float(x) for x in line.split()[:3]])

    return DFTResult(
        energy=energy,
        forces=np.array(forces),
        stress=np.array(stress)
    )
