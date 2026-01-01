# src/mlip_autopipec/utils/dft_utils.py

import re
from typing import Any

import numpy as np


def parse_qe_output(output: str) -> dict[str, Any]:
    """
    Parses the output of a Quantum Espresso calculation to extract key results.

    Args:
        output: The stdout content from the Quantum Espresso run.

    Returns:
        A dictionary containing the final energy, forces, and stress tensor.
        Returns an empty dictionary if parsing fails.
    """
    results = {}

    # Find the final energy
    energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
    if energy_match:
        # Convert energy from Rydberg to eV
        results["energy"] = float(energy_match.group(1)) * 13.605693122994

    # Find the forces using a block-based approach for robustness
    force_block_match = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):"
        r"(.*?)(?:Total force|total\s+stress)",
        output,
        flags=re.DOTALL,  # Use re.DOTALL to make '.' match newlines
    )
    if force_block_match:
        force_block_text = force_block_match.group(1)
        forces = []
        for line in force_block_text.strip().split("\n"):
            line = line.strip()
            if line.startswith("atom"):
                parts = line.split()
                if len(parts) >= 7:
                    forces.append([float(parts[4]), float(parts[5]), float(parts[6])])
        if forces:
            # Convert forces from Ry/au to eV/Angstrom
            ry_au_to_ev_a = 13.605693122994 / 0.529177210903
            results["forces"] = np.array(forces) * ry_au_to_ev_a

    # Find the stress tensor
    stress_match = re.search(
        r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*(-?\d+\.\d+)\s*\n"
        r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n)"
        r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n)"
        r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n)",
        output,
    )
    if stress_match:
        lines = stress_match.group(0).strip().split("\n")[1:]
        stress = np.array([list(map(float, line.split())) for line in lines])
        # Convert stress from kbar to eV/Angstrom^3
        kbar_to_ev_a3 = 1 / 160.21766208
        results["stress"] = stress * kbar_to_ev_a3

    return results
