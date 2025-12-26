import re
from typing import Dict, Any, Tuple, Optional
import numpy as np
from ase.atoms import Atoms

def generate_qe_input(atoms: Atoms, parameters: Dict[str, Any]) -> str:
    """Generates a Quantum Espresso input file string for a given ASE Atoms object."""

    input_lines = []

    # Control block
    input_lines.append("&CONTROL")
    input_lines.append("    calculation = 'scf'")
    input_lines.append("    restart_mode = 'from_scratch'")
    input_lines.append(f"    pseudo_dir = '{parameters.get('pseudo_dir', '.')}'")
    input_lines.append("    outdir = './out'")
    input_lines.append("    prefix = 'qe_run'")
    input_lines.append("    tprnfor = .true.")
    input_lines.append("    tstress = .true.")
    input_lines.append("/")
    input_lines.append("")

    # System block
    input_lines.append("&SYSTEM")
    input_lines.append(f"    ibrav = 0")
    input_lines.append(f"    nat = {len(atoms)}")
    input_lines.append(f"    ntyp = {len(set(atoms.get_chemical_symbols()))}")
    input_lines.append(f"    ecutwfc = {parameters.get('ecutwfc', 60.0)}")
    input_lines.append("/")
    input_lines.append("")

    # Electrons block
    input_lines.append("&ELECTRONS")
    input_lines.append("    mixing_beta = 0.7")
    input_lines.append("    conv_thr = 1.0e-8")
    input_lines.append("/")
    input_lines.append("")

    # Atomic species
    input_lines.append("ATOMIC_SPECIES")
    species = sorted(list(set(atoms.get_chemical_symbols())))
    mass_map = {s: atoms[atoms.get_chemical_symbols().index(s)].mass for s in species}
    for symbol in species:
        pseudo = parameters['pseudopotentials'][symbol]
        input_lines.append(f"    {symbol} {mass_map[symbol]:.4f} {pseudo}")
    input_lines.append("")

    # Atomic positions
    input_lines.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        input_lines.append(f"    {atom.symbol} {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}")
    input_lines.append("")

    # K-points
    input_lines.append("K_POINTS {automatic}")
    kpoints = parameters.get('kpoints', [1, 1, 1, 0, 0, 0])
    input_lines.append(f"    {' '.join(map(str, kpoints))}")
    input_lines.append("")

    # Cell parameters
    input_lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.cell:
        input_lines.append(f"    {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")

    return "\n".join(input_lines)


def parse_qe_output(output: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Parses Quantum Espresso output to extract key information."""

    result = {}
    error_message = None

    if "JOB DONE" not in output:
        return None, "Calculation did not finish."

    energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
    if energy_match:
        result['total_energy_ev'] = float(energy_match.group(1)) * 13.6057
    else:
        return None, "Total energy not found."

    forces_match = re.search(r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n(.*?)\n\s*\n", output, re.DOTALL)
    if forces_match:
        forces_block = forces_match.group(1)
        force_lines = re.findall(r"atom\s+\d+\s+type\s+\w+\s+force\s+=\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", forces_block)
        forces = np.array(force_lines, dtype=float) * 25.71104
        result['forces'] = forces.tolist()
    else:
        return None, "Forces block not found."

    stress_match = re.search(r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+pbar\s*\n((?:\s*[-.\d]+\s*){9})", output)
    if stress_match:
        stress_values = np.fromstring(stress_match.group(1), sep=' ')
        stress_kbar = stress_values.reshape(3, 3)
        result['stress'] = (stress_kbar * 0.1).tolist() # Convert kbar to GPa
    else:
        return None, "Stress block not found."

    return result, error_message
