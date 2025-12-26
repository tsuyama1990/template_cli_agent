import re
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTResult

# Conversion factors
RY_TO_EV = 13.605693122994
RY_AU_TO_EV_A = 25.711043 * 2 # Rydbergs per Bohr to eV per Angstrom

def generate_qe_input(atoms: Atoms) -> str:
    """
    Generates a basic Quantum Espresso input file for an SCF calculation.

    Args:
        atoms: The atomic structure.

    Returns:
        A string containing the QE input file content.
    """
    species = sorted(list(set(atoms.get_chemical_symbols())))
    pseudopotentials = {spec: f"{spec}.upf" for spec in species}

    input_lines = [
        "&CONTROL",
        "    calculation = 'scf'",
        "    restart_mode = 'from_scratch'",
        "    prefix = 'pwscf'",
        "    outdir = './out'",
        "    pseudo_dir = './pseudos'",
        "/",
        "&SYSTEM",
        f"    ibrav = 0",
        f"    nat = {len(atoms)}",
        f"    ntyp = {len(species)}",
        "    ecutwfc = 60.0",
        "/",
        "&ELECTRONS",
        "    mixing_beta = 0.7",
        "    conv_thr = 1.0e-8",
        "/",
        "ATOMIC_SPECIES",
    ]
    for spec in species:
        # A placeholder mass is used, as QE uses it from the pseudo.
        input_lines.append(f"    {spec}  28.0855  {pseudopotentials[spec]}")

    input_lines.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        input_lines.append(f"    {atom.symbol}  {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}")

    input_lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.cell:
        input_lines.append(f"    {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")

    input_lines.append("K_POINTS {automatic}")
    input_lines.append("    6 6 6 0 0 0") # A reasonable default

    return "\n".join(input_lines)


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the output of a Quantum Espresso calculation to extract key results.

    Args:
        output: The stdout from the 'pw.x' run.

    Returns:
        A DFTResult object containing the parsed data.
    """
    total_energy_match = re.search(r"!\s+total energy\s+=\s+([-.\d]+)\s+Ry", output)
    forces_match = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+([-.\d\s]+)", output)
    stress_match = re.search(r"total stress\s+\(Ry/bohr\*\*3\)\s+[\s\S]*?P=\s*[-.\d]+\n((?:[-.\d\s]+\n){3})", output)

    # Check for common error messages
    error_patterns = [
        r"S matrix not positive definite",
        r"SCF NOT CONVERGED",
        r"convergence has not been achieved"
    ]
    for pattern in error_patterns:
        if re.search(pattern, output, re.IGNORECASE):
            return DFTResult(
                total_energy_ev=0.0,
                forces=[],
                stress=[],
                was_successful=False,
                error_message=re.search(pattern, output, re.IGNORECASE).group(0),
            )

    if not total_energy_match:
         return DFTResult(
                total_energy_ev=0.0,
                forces=[],
                stress=[],
                was_successful=False,
                error_message="Could not find total energy in QE output.",
            )

    total_energy = float(total_energy_match.group(1)) * RY_TO_EV

    forces = []
    if forces_match:
        for force_line in forces_match:
            forces.append([float(f) * RY_AU_TO_EV_A for f in force_line.split()])

    stress = []
    if stress_match:
         stress_lines = stress_match.group(1).strip().split('\n')
         for line in stress_lines:
             stress.append([float(s) for s in line.split()])


    return DFTResult(
        total_energy_ev=total_energy,
        forces=forces,
        stress=stress,
        was_successful=True,
    )
