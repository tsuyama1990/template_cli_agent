import re

from ase.atoms import Atoms


def generate_qe_input(atoms: Atoms, parameters: dict) -> str:
    """
    Generates a Quantum Espresso input file string from an ASE Atoms object.
    """
    pseudos = {atom.symbol: f"{atom.symbol}.UPF" for atom in atoms}
    ecutwfc = parameters.get("ecutwfc", 60)
    ecutrho = parameters.get("ecutrho", 240)
    k_points = parameters.get("k_points", [1, 1, 1])
    pseudo_dir = parameters.get("pseudo_dir", "./pseudos")


    input_str = f"""
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'pwscf'
    outdir = './out'
    pseudo_dir = '{pseudo_dir}'
/
&SYSTEM
    ibrav = 0
    nat = {len(atoms)}
    ntyp = {len(set(atoms.get_chemical_symbols()))}
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
/
&ELECTRONS
    mixing_beta = 0.7
    conv_thr = 1.0e-8
/
ATOMIC_SPECIES
"""
    for symbol in sorted(set(atoms.get_chemical_symbols())):
        input_str += f"  {symbol} 1.0 {pseudos[symbol]}\n"

    input_str += "\nATOMIC_POSITIONS {angstrom}\n"
    for atom in atoms:
        pos = atom.position
        input_str += f"  {atom.symbol} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n"

    input_str += "\nCELL_PARAMETERS {angstrom}\n"
    for vector in atoms.get_cell():
        input_str += f"  {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}\n"

    input_str += "\nK_POINTS {automatic}\n"
    input_str += f"  {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0\n"

    return input_str


def parse_qe_output(
    output: str,
) -> tuple[float | None, list[list[float]] | None, list[list[float]] | None, bool, str | None]:
    """
    Parses the output of a Quantum Espresso calculation to extract energy, forces, and stress.
    """
    energy_match = re.search(r"!    total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
    if not energy_match:
        return None, None, None, False, "Could not find total energy."

    ry_to_ev = 13.605693122994
    energy = float(energy_match.group(1)) * ry_to_ev

    force_block_match = re.search(
        r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n(.+?)\n\n", output, re.DOTALL
    )
    if not force_block_match:
        return energy, None, None, False, "Could not find forces."

    ry_au_to_ev_a = ry_to_ev / 0.529177210903
    forces = []
    force_lines = force_block_match.group(1).strip().split('\n')
    for line in force_lines:
        parts = line.split()
        fx, fy, fz = map(float, parts[-3:])
        forces.append([fx * ry_au_to_ev_a, fy * ry_au_to_ev_a, fz * ry_au_to_ev_a])

    stress_block_match = re.search(
        r"total stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*(-?\d+\.\d+)\n((?:\s+-?\d+\.\d+){6})",
        output,
    )
    if not stress_block_match:
        return energy, forces, None, False, "Could not find stress."

    kbar_to_ev_a3 = 1.0 / 160.21766208
    stress_values = list(map(float, stress_block_match.group(2).strip().split()))
    stress_matrix = [
        [
            stress_values[0] * kbar_to_ev_a3,
            stress_values[3] * kbar_to_ev_a3,
            stress_values[4] * kbar_to_ev_a3,
        ],
        [
            stress_values[3] * kbar_to_ev_a3,
            stress_values[1] * kbar_to_ev_a3,
            stress_values[5] * kbar_to_ev_a3,
        ],
        [
            stress_values[4] * kbar_to_ev_a3,
            stress_values[5] * kbar_to_ev_a3,
            stress_values[2] * kbar_to_ev_a3,
        ],
    ]

    return energy, forces, stress_matrix, True, None
