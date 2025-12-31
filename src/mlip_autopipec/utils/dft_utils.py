import re
import numpy as np
from ase.stress import full_3x3_to_voigt_6_stress


def parse_qe_output(output: str) -> dict:
    """
    Parses the output of a Quantum Espresso calculation to extract energy,
    forces, and stress.

    Args:
        output: The stdout content from the QE calculation.

    Returns:
        A dictionary containing the extracted data.
    """
    data = {}

    # More robust regex for forces, looking for the "Total force" line as a delimiter.
    force_match = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry\/au\):\s*\n(.*?)(?=\n\s*The total force is)",
        output,
        re.DOTALL,
    )
    if force_match:
        force_lines = force_match.group(1).strip().split("\n")
        forces = [list(map(float, line.split()[-3:])) for line in force_lines]
        # Convert from Ry/au to eV/Angstrom
        data["forces"] = np.array(forces) * (13.605693122994 / 0.529177210903)

    stress_match = re.search(
        r"total stress\s*\(Ry\/bohr\*\*3\)\s*\n(.*?)(?=\n\n|Writing output data file)",
        output,
        re.DOTALL,
    )
    if stress_match:
        stress_lines = stress_match.group(1).strip().split("\n")
        # Ensure we only parse the 3x3 matrix part
        stress_matrix = np.array(
            [list(map(float, line.split()[:3])) for line in stress_lines]
        )

        if stress_matrix.shape == (3, 3):
            # Convert from Ry/bohr^3 to eV/Ang^3
            ry_to_ev = 13.605693122994
            bohr_to_ang = 0.529177210903
            conversion_factor = ry_to_ev / (bohr_to_ang**3)
            stress_matrix_ev_ang3 = stress_matrix * conversion_factor

            # Convert the full 3x3 tensor to Voigt form for ASE
            data["stress"] = full_3x3_to_voigt_6_stress(stress_matrix_ev_ang3)

    energy_match = re.search(r"!\s*total energy\s*=\s*([-\d.]+)\s*Ry", output)
    if energy_match:
        # Convert from Ry to eV
        data["energy"] = float(energy_match.group(1)) * 13.605693122994

    return data
