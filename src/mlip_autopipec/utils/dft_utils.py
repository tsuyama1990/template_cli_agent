import re
import numpy as np

# Conversion factors from ASE
RY_TO_EV = 13.605693122994
RY_AU_TO_EV_ANG = RY_TO_EV / 0.529177210903
KBAR_TO_EV_ANG3 = 1.0 / 160.21766208


def parse_qe_output(output_string: str) -> dict:
    """
    Parses the output of a Quantum Espresso pw.x calculation to extract
    final energy, forces, and stress.

    Args:
        output_string: The stdout of the pw.x calculation as a string.

    Returns:
        A dictionary containing the parsed energy, forces, and stress,
        converted to ASE's standard units (eV, eV/Ang, eV/Ang^3).
    """
    # --- Energy ---
    energy_match = re.search(r"!\s+total energy\s+=\s+([-\d\.]+)\s+Ry", output_string)
    if not energy_match:
        raise ValueError("Could not find total energy in QE output.")
    energy = float(energy_match.group(1)) * RY_TO_EV

    # --- Forces ---
    # CORRECTED: Added \s* to handle leading whitespace before "Total Force"
    force_block_match = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n(.+?)\n\n\s*Total Force",
        output_string,
        re.DOTALL,
    )
    if not force_block_match:
        raise ValueError("Could not find forces block in QE output.")

    forces = []
    force_lines = force_block_match.group(1).strip().split("\n")
    for line in force_lines:
        parts = line.split()
        fx, fy, fz = map(float, parts[-3:])
        forces.append([fx, fy, fz])
    forces_array = np.array(forces) * RY_AU_TO_EV_ANG

    # --- Stress ---
    # The stress block in QE is in kbar, but ASE wants eV/Ang^3.
    # The header can be "total stress" or "The total stress is."
    # Let's find the header first, then the block.
    stress_header_match = re.search(
        r"(?:total stress|The total stress is)\.\s*\(Ry/bohr\*\*3\)",
        output_string
    )
    if not stress_header_match:
        raise ValueError("Could not find stress block header in QE output.")

    # Find the 3x3 matrix following the header
    stress_block_match = re.search(
        r"(-?[\d\.]+\s+-?[\d\.]+\s+-?[\d\.]+\s*\n\s*-?[\d\.]+\s+-?[\d\.]+\s+-?[\d\.]+\s*\n\s*-?[\d\.]+\s+-?[\d\.]+\s+-?[\d\.]+)",
        output_string[stress_header_match.end():]
    )
    if not stress_block_match:
        raise ValueError("Could not find stress block matrix in QE output.")

    stress_flat = [float(s) for s in stress_block_match.group(1).split()]
    stress_voigt = np.array(stress_flat)

    # QE prints in Voigt order: (xx, yy, zz, xy, xz, yz)
    # but the matrix is formatted differently, let's read it line by line
    stress_lines = stress_block_match.group(1).strip().split('\n')
    stress_matrix = np.array([list(map(float, line.split())) for line in stress_lines])

    # Convert from Ry/bohr^3 to eV/Ang^3 (ASE's internal unit)
    stress_matrix_ev_ang3 = stress_matrix * (RY_TO_EV / (0.529177**3))


    return {
        "energy": energy,
        "forces": forces_array,
        "stress": stress_matrix_ev_ang3,
    }
