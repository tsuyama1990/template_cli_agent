"""Utilities for interacting with Quantum Espresso."""
import re
import numpy as np
from typing import Dict, Any
from ase.atoms import Atoms
from ..data.models import DFTResult

# Conversion factors
RY_TO_EV = 13.605693122994
AU_TO_ANGSTROM = 0.529177210903
RY_AU_TO_EV_ANG = RY_TO_EV / AU_TO_ANGSTROM
KBAR_TO_EV_ANG3 = 1 / 1602.1766208

def _format_qe_value(value: Any) -> str:
    """Formats a Python value into a QE-compatible string."""
    if isinstance(value, bool):
        return f".{str(value).lower()}."
    if isinstance(value, str) and value.startswith("'") and value.endswith("'"):
        return value
    if isinstance(value, str):
        return f"'{value}'"
    return str(value)

def generate_qe_input(atoms: Atoms, parameters: Dict[str, Any], pseudos: Dict[str, str]) -> str:
    """
    Generates a Quantum Espresso input file content from an ASE Atoms object.

    Args:
        atoms: The ASE Atoms object.
        parameters: Dictionary of QE input parameters.
        pseudos: Dictionary mapping atom species to pseudopotential filenames.

    Returns:
        A string containing the QE input file content.
    """

    input_parts = []

    # &CONTROL, &SYSTEM, &ELECTRONS
    for namelist in ["CONTROL", "SYSTEM", "ELECTRONS"]:
        params = parameters.get(namelist.lower(), {})
        input_parts.append(f"&{namelist}")
        if namelist == "SYSTEM":
            input_parts.append(f"  ibrav = 0")
            input_parts.append(f"  nat = {len(atoms)}")
            input_parts.append(f"  ntyp = {len(set(atoms.get_chemical_symbols()))}")

        for key, value in params.items():
            input_parts.append(f"  {key} = {_format_qe_value(value)}")
        input_parts.append("/")

    # ATOMIC_SPECIES
    input_parts.append("ATOMIC_SPECIES")
    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    for symbol in unique_symbols:
        input_parts.append(f"  {symbol}  1.0  {pseudos[symbol]}")

    # ATOMIC_POSITIONS
    input_parts.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        input_parts.append(f"  {atom.symbol}  {atom.position[0]:.8f}  {atom.position[1]:.8f}  {atom.position[2]:.8f}")

    # K_POINTS
    k_points = parameters.get("k_points", {"scheme": "automatic", "grid": [1, 1, 1]})
    input_parts.append(f"K_POINTS {{{k_points['scheme']}}}")
    if k_points['scheme'] == 'automatic':
        grid = k_points['grid']
        shifts = k_points.get('shifts', [0, 0, 0])
        input_parts.append(f"  {grid[0]} {grid[1]} {grid[2]} {shifts[0]} {shifts[1]} {shifts[2]}")

    # CELL_PARAMETERS
    input_parts.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.cell:
        input_parts.append(f"  {vector[0]:.8f}  {vector[1]:.8f}  {vector[2]:.8f}")

    return "\n".join(input_parts) + "\n"


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the output of a Quantum Espresso calculation to extract results.

    Args:
        output: The stdout content from the 'pw.x' run.

    Returns:
        A DFTResult object containing the parsed data.
    """
    was_successful = "JOB DONE" in output
    error_message = None
    if not was_successful:
        if "convergence has not been achieved" in output:
            error_message = "SCF did not converge"
        else:
            error_match = re.search(r"Error in routine .*?\n\s*(.*?)\n", output)
            if error_match:
                error_message = error_match.group(1).strip()
            else:
                error_message = "Unknown error: 'JOB DONE' not found."

    energy_matches = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
    total_energy_ev = float(energy_matches[-1]) * RY_TO_EV if energy_matches else 0.0

    lines = output.split('\n')
    num_atoms_match = re.search(r'number of atoms/cell\s+=\s+(\d+)', output)
    num_atoms = int(num_atoms_match.group(1)) if num_atoms_match else 0

    forces = [[0.0] * 3] * num_atoms
    stress = [[0.0] * 3] * 3


    # Find and parse forces
    try:
        force_section_found = False
        for i, line in enumerate(lines):
            if "Forces acting on atoms" in line:
                force_section_found = True
                for force_line in lines[i+1:]:
                    force_line_stripped = force_line.strip()
                    if not force_line_stripped:
                        continue
                    if "atom" not in force_line_stripped:
                        break
                    parts = force_line_stripped.split()
                    fx, fy, fz = map(float, parts[6:9])
                    forces.append([fx * RY_AU_TO_EV_ANG, fy * RY_AU_TO_EV_ANG, fz * RY_AU_TO_EV_ANG])
                break
        if not force_section_found and was_successful and num_atoms > 0:
             forces = [[0.0, 0.0, 0.0]] * num_atoms

    except (ValueError, IndexError):
        if was_successful:
            was_successful = False
            error_message = "Failed to parse forces from a supposedly successful QE run."
        forces = []

    # Find and parse stress
    try:
        stress_section_found = False
        for i, line in enumerate(lines):
            if "total stress" in line and "(Ry/bohr**3)" in line:
                stress_section_found = True
                stress_lines = lines[i+1:i+4]
                stress_matrix_ry_bohr3 = [list(map(float, l.split()[:3])) for l in stress_lines]
                stress = (np.array(stress_matrix_ry_bohr3) * RY_TO_EV / (AU_TO_ANGSTROM**3)).tolist()
                break
    except (ValueError, IndexError):
        if was_successful:
            was_successful = False
            error_message = "Failed to parse stress from a supposedly successful QE run."
        stress = [[0.0] * 3] * 3

    return DFTResult(
        total_energy_ev=total_energy_ev,
        forces=forces,
        stress=stress,
        was_successful=was_successful,
        error_message=error_message,
    )
