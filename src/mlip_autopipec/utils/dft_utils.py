import re
import numpy as np
from ase.atoms import Atoms
from typing import Dict, Any, Tuple

from mlip_autopipec.data.models import DFTResult

# Conversion factors
RY_TO_EV = 13.605693122994
BOHR_TO_ANGSTROM = 0.529177210903
RY_AU_TO_EV_A = RY_TO_EV / BOHR_TO_ANGSTROM

def generate_qe_input(atoms: Atoms, params: Dict[str, Any], pseudos: Dict[str, str]) -> str:
    """Generates a Quantum Espresso input file string."""
    nat = len(atoms)
    ntyp = len(pseudos)

    control = params.get('control', {})
    system = params.get('system', {})
    electrons = params.get('electrons', {})

    input_lines = [
        "&CONTROL",
        f"    calculation = '{control.get('calculation', 'scf')}'",
        f"    restart_mode = '{control.get('restart_mode', 'from_scratch')}'",
        f"    prefix = '{control.get('prefix', 'pwscf')}'",
        "    tstress = .true.",
        "    tprnfor = .true.",
        f"    outdir = '{control.get('outdir', './')}'",
        f"    pseudo_dir = '{control.get('pseudo_dir', './')}'",
        "/",
        "&SYSTEM",
        f"    ibrav = 0",
        f"    nat = {nat}",
        f"    ntyp = {ntyp}",
        f"    ecutwfc = {system.get('ecutwfc', 60.0)}",
        f"    ecutrho = {system.get('ecutrho', 240.0)}",
        "/",
        "&ELECTRONS",
        f"    mixing_beta = {electrons.get('mixing_beta', 0.7)}",
        f"    conv_thr = {electrons.get('conv_thr', 1.0e-10)}",
        "/",
        "ATOMIC_SPECIES",
    ]

    unique_symbols = sorted(list(pseudos.keys()))
    for symbol in unique_symbols:
        input_lines.append(f"    {symbol} {atoms.get_masses()[list(atoms.get_chemical_symbols()).index(symbol)]} {pseudos[symbol]}")

    input_lines.append("\nCELL_PARAMETERS angstrom")
    for vector in atoms.get_cell():
        input_lines.append(f"    {vector[0]:.10f} {vector[1]:.10f} {vector[2]:.10f}")

    input_lines.append("\nATOMIC_POSITIONS angstrom")
    for atom in atoms:
        input_lines.append(f"    {atom.symbol} {atom.position[0]:.10f} {atom.position[1]:.10f} {atom.position[2]:.10f}")

    k_points = params.get('k_points', {'scheme': 'automatic', 'grid': [1, 1, 1], 'shift': [0, 0, 0]})
    input_lines.append(f"\nK_POINTS {k_points['scheme']}")
    if k_points['scheme'] == 'automatic':
        grid = k_points['grid']
        shift = k_points['shift']
        input_lines.append(f"  {grid[0]} {grid[1]} {grid[2]} {shift[0]} {shift[1]} {shift[2]}")

    return "\n".join(input_lines)


def parse_qe_output(output: str) -> DFTResult:
    """Parses Quantum Espresso output to extract energy, forces, and stress."""
    total_energy_ev = None
    forces = []
    stress = np.zeros((3, 3))
    was_successful = "JOB DONE" in output
    error_message = None

    if not was_successful:
        if "SCF NOT CONVERGED" in output:
            error_message = "SCF did not converge."
        else:
            error_message = "QE job did not finish successfully."

    # Find total energy
    energy_match = re.search(r"!    total energy\s+=\s+(-?[\d\.]+)\s+Ry", output)
    if energy_match:
        total_energy_ev = float(energy_match.group(1)) * RY_TO_EV

    # Find forces
    force_match = re.search(r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n\s*\n((?:\s*atom\s+\d+\s+type\s+\w+\s+force\s+=\s+[-]?[\d\.]+\s+[-]?[\d\.]+\s+[-]?[\d\.]+\s*\n)+)", output)
    if force_match:
        force_lines = force_match.group(1).strip().split('\n')
        for line in force_lines:
            parts = line.split('=')
            force_values = [float(f) for f in parts[1].strip().split()]
            forces.append([f * RY_AU_TO_EV_A for f in force_values])

    # Find stress
    stress_match = re.search(r"total stress\s*\(Ry/bohr\*\*3\)\s*\n((?:\s*[-]?[\d\.]+\s*[-]?[\d\.]+\s*[-]?[\d\.]+\s*\n){3})", output)
    if stress_match:
        stress_lines = stress_match.group(1).strip().split('\n')
        stress_ry_bohr3 = np.array([list(map(float, line.split())) for line in stress_lines])
        # Convert Ry/bohr^3 to eV/Angstrom^3
        stress = stress_ry_bohr3 * RY_TO_EV / (BOHR_TO_ANGSTROM**3)

    if total_energy_ev is None or not forces:
        was_successful = False
        if error_message is None:
            error_message = "Could not parse energy or forces from output."

    return DFTResult(
        total_energy_ev=total_energy_ev if total_energy_ev is not None else 0.0,
        forces=forces if forces else [],
        stress=stress.tolist() if stress is not None else [],
        was_successful=was_successful,
        error_message=error_message
    )
