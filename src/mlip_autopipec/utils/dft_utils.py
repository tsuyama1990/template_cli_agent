import re
from typing import Any, Dict

import numpy as np
from ase.atoms import Atoms


def generate_qe_input(
    atoms: Atoms, parameters: Dict[str, Any], pseudopotentials: Dict[str, str]
) -> str:
    """
    Generates a Quantum Espresso input file content.

    Args:
        atoms: The ASE Atoms object.
        parameters: A dictionary of QE parameters.
        pseudopotentials: A mapping from atomic symbol to pseudopotential filename.

    Returns:
        The content of the QE input file as a string.
    """
    input_lines = []

    # &CONTROL
    input_lines.append("&CONTROL")
    control_params = parameters.get("control", {})
    control_params.setdefault("calculation", "scf")
    control_params.setdefault("verbosity", "high")
    for key, value in control_params.items():
        if isinstance(value, str):
            input_lines.append(f"  {key} = '{value}'")
        else:
            input_lines.append(f"  {key} = {value}")
    input_lines.append("/")

    # &SYSTEM
    input_lines.append("&SYSTEM")
    system_params = parameters.get("system", {})
    system_params.setdefault("ibrav", 0)
    system_params.setdefault("nat", len(atoms))
    system_params.setdefault("ntyp", len(set(atoms.get_chemical_symbols())))
    for key, value in system_params.items():
        input_lines.append(f"  {key} = {value}")
    input_lines.append("/")

    # &ELECTRONS
    input_lines.append("&ELECTRONS")
    electron_params = parameters.get("electrons", {})
    for key, value in electron_params.items():
        input_lines.append(f"  {key} = {value}")
    input_lines.append("/")

    # ATOMIC_SPECIES
    input_lines.append("ATOMIC_SPECIES")
    species = sorted(list(set(atoms.get_chemical_symbols())))
    masses = {symbol: atoms[i].mass for i, symbol in enumerate(atoms.get_chemical_symbols())}
    for symbol in species:
        input_lines.append(f"  {symbol} {masses[symbol]:.4f} {pseudopotentials[symbol]}")

    # ATOMIC_POSITIONS
    input_lines.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        pos = atom.position
        input_lines.append(f"  {atom.symbol} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}")

    # K_POINTS
    input_lines.append("K_POINTS {automatic}")
    k_points = parameters.get("k_points", [1, 1, 1, 0, 0, 0])
    input_lines.append(f"  {' '.join(map(str, k_points))}")

    # CELL_PARAMETERS
    input_lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.get_cell():
        input_lines.append(f"  {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")

    return "\n".join(input_lines)


def parse_qe_output(output: str) -> tuple[bool, str | None, dict[str, Any] | None]:
    """
    Parses the output of a Quantum Espresso calculation.

    Args:
        output: The stdout content from the pw.x run.

    Returns:
        A tuple containing:
        - was_successful (bool): True if the run converged and results were found.
        - error_message (str | None): An error message if not successful.
        - results (dict | None): A dictionary with 'energy', 'forces', and 'stress' if successful.
    """
    if "JOB DONE" not in output:
        return False, "Calculation did not finish (JOB DONE not found).", None
    if "convergence has been achieved" not in output and "convergence NOT achieved" in output:
        return False, "SCF calculation did not converge.", None

    results = {}
    try:
        # Energy (convert from Ry to eV)
        energy_ry = float(re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output).group(1))
        results["total_energy_ev"] = energy_ry * 13.605693122994

        # Forces (convert from Ry/au to eV/Angstrom)
        force_section = re.search(
            r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n([\s\S]+?)\n\n", output
        )
        forces_ry_au = []
        for line in force_section.group(1).splitlines():
            parts = line.split()
            forces_ry_au.append([float(parts[6]), float(parts[7]), float(parts[8])])
        results["forces"] = (np.array(forces_ry_au) * (13.605693122994 / 0.529177210903)).tolist()

        # Stress (convert from kbar to eV/Angstrom^3)
        stress_section = re.search(
            r"total stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*[\d\.\-]+\n([\s\S]+?)(?:\n\n|!|$)",
            output,
        )
        stress_kbar_lines = stress_section.group(1).strip().splitlines()
        # The stress block in QE output has 6 columns, but we only need the first 3 for the tensor
        stress_kbar = np.array([list(map(float, line.split()[:3])) for line in stress_kbar_lines])
        # 1 kbar = 1/160.21766208 GPa; 1 eV/A^3 = 160.21766208 GPa
        results["stress"] = (stress_kbar / 160.21766208).tolist()

    except (AttributeError, IndexError, ValueError) as e:
        return False, f"Failed to parse output: {e}", None

    return True, None, results
