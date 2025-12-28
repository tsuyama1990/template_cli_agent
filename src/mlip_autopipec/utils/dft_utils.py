# Description: Utility functions for interacting with Quantum Espresso.
import re

import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTResult

# Basic SSSP-like recommendations for a few elements for demonstration
PSEUDO_POTENTIALS = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
ECUTWFC = {"Si": 60}


def generate_qe_input(atoms: Atoms, parameters: dict | None = None) -> str:
    """
    Generates a Quantum Espresso input file string for a given ASE Atoms object.

    Args:
        atoms: The ASE Atoms object for the structure.
        parameters: A dictionary of QE parameters to override defaults.

    Returns:
        A string containing the content of the QE input file.
    """
    if parameters is None:
        parameters = {}

    # System parameters
    nat = len(atoms)
    ntyp = len(set(atoms.get_chemical_symbols()))
    ibrav = 0  # Use explicit lattice vectors

    # Default parameters, can be overridden
    control_params = {
        "calculation": "'scf'",
        "verbosity": "'high'",
        "prefix": "'mlip_autopipec'",
        "pseudo_dir": "'./'",
        "outdir": "'./out/'",
    }
    system_params = {
        "nat": nat,
        "ntyp": ntyp,
        "ibrav": ibrav,
        "ecutwfc": ECUTWFC.get(atoms.get_chemical_symbols()[0], 30),
    }
    electrons_params = {
        "mixing_beta": 0.7,
        "conv_thr": 1.0e-10,
    }

    # Update with any user-provided parameters
    control_params.update(parameters.get("CONTROL", {}))
    system_params.update(parameters.get("SYSTEM", {}))
    electrons_params.update(parameters.get("ELECTRONS", {}))

    # Build the input string
    lines = []
    lines.append("&CONTROL")
    for key, value in control_params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")

    lines.append("&SYSTEM")
    for key, value in system_params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")

    lines.append("&ELECTRONS")
    for key, value in electrons_params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")

    # Atomic species
    lines.append("ATOMIC_SPECIES")
    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    for symbol in unique_symbols:
        mass = atoms[atoms.symbols.index(symbol)].mass
        pseudo = PSEUDO_POTENTIALS.get(symbol, f"{symbol}.UPF_placeholder")
        lines.append(f"  {symbol}   {mass:.4f}   {pseudo}")

    # Atomic positions
    lines.append("ATOMIC_POSITIONS {crystal}")
    for atom in atoms:
        lines.append(
            f"  {atom.symbol} {atom.scaled_position[0]:.8f} "
            f"{atom.scaled_position[1]:.8f} {atom.scaled_position[2]:.8f}"
        )

    # K-points
    lines.append("K_POINTS {automatic}")
    kpts = parameters.get("K_POINTS", [4, 4, 4, 0, 0, 0])
    lines.append(f"  {' '.join(map(str, kpts))}")

    # Cell parameters
    lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.cell:
        lines.append(f"  {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")

    return "\n".join(lines)


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the output of a Quantum Espresso calculation to extract results.

    Args:
        output: A string containing the stdout of the QE run.

    Returns:
        A DFTResult object containing the parsed data.
    """
    total_energy_ev: float | None = None
    forces: list[list[float]] = []
    stress: list[list[float]] = []
    was_successful = False
    error_message: str | None = None

    try:
        # Check for successful completion
        if "JOB DONE." in output:
            was_successful = True

        # Check for convergence issues
        if "convergence NOT achieved" in output:
            was_successful = False
            error_message = "SCF did not converge."

        # Parse total energy
        energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
        if energy_match:
            ry_to_ev = 13.605693122994
            total_energy_ev = float(energy_match.group(1)) * ry_to_ev

        # Parse forces
        force_section_match = re.search(
            r"Forces acting on atoms\s*\(Ry/au\):\s*\n(.*?)(?=\n\s*\n|\n\s*Total force|\n\s*total stress|\n\s*! total energy|\n\s*JOB DONE\.)",
            output,
            re.DOTALL,
        )
        if force_section_match:
            force_lines = force_section_match.group(1).strip().split("\n")
            for line in force_lines:
                parts = line.split()
                if len(parts) >= 9 and parts[0] == "atom":
                    force_ry_au = [float(parts[6]), float(parts[7]), float(parts[8])]
                    # Convert Ry/a.u. to eV/Angstrom
                    ry_au_to_ev_ang = 25.71104309541616 * 2
                    force_ev_ang = [f * ry_au_to_ev_ang for f in force_ry_au]
                    forces.append(force_ev_ang)

        # Parse stress tensor
        stress_section_match = re.search(
            r"total stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*(-?\d+\.\d+)\n"
            r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\n"
            r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\n"
            r"\s*(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)",
            output,
        )
        if stress_section_match:
            ry_bohr3_to_ev_ang3 = 13.605693122994 / (0.529177210903 ** 3)
            stress_lines = [
                stress_section_match.group(i) for i in range(2, 5)
            ]  # groups 2, 3, 4
            stress_tensor_ry = [list(map(float, line.split())) for line in stress_lines]
            stress = (
                (np.array(stress_tensor_ry) * ry_bohr3_to_ev_ang3).tolist()
                if stress_tensor_ry
                else []
            )

        # If successful but key data is missing, mark as failed
        if was_successful and (total_energy_ev is None or not forces):
            was_successful = False
            error_message = "Parsing failed: Could not find energy or forces in output."

    except Exception as e:
        was_successful = False
        error_message = f"An unexpected parsing error occurred: {e}"

    # Generic error if parsing fails without a specific reason
    if not was_successful and error_message is None:
        error_message = "Malformed or empty QE output."

    return DFTResult(
        total_energy_ev=total_energy_ev if total_energy_ev is not None else 0.0,
        forces=forces,
        stress=stress,
        was_successful=was_successful,
        error_message=error_message,
    )
