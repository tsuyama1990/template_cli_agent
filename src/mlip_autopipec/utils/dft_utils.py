import re
from typing import Any

from ase import Atoms
from ase.data import atomic_masses

RY_TO_EV = 13.605693122994
KBAR_TO_EV_PER_ANG3 = 1 / 160.21766208
RY_BOHR3_TO_EV_ANG3 = RY_TO_EV / (0.529177**3)

def generate_qe_input(
    atoms: Atoms,
    parameters: dict[str, Any],
    pseudopotentials: dict[str, str],
    kpoints: tuple[int, int, int],
) -> str:
    """
    Generates a Quantum Espresso input file string from an ASE Atoms object and parameters.
    """
    # &CONTROL block
    control_block = (
        "&CONTROL\n"
        f"    calculation = '{parameters.get('calculation', 'scf')}'\n"
        "    restart_mode = 'from_scratch'\n"
        "    prefix = 'pwscf'\n"
        "    outdir = './out'\n"
        "    pseudo_dir = './pseudos'\n"
        "/\n"
    )

    # &SYSTEM block
    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    ntyp = len(unique_symbols)
    system_block = (
        "&SYSTEM\n"
        "    ibrav = 0\n"
        f"    nat = {len(atoms)}\n"
        f"    ntyp = {ntyp}\n"
        f"    ecutwfc = {parameters['ecutwfc']}\n"
        f"    ecutrho = {parameters['ecutrho']}\n"
        f"    occupations = '{parameters.get('occupations', 'smearing')}'\n"
        f"    smearing = '{parameters.get('smearing', 'gaussian')}'\n"
        f"    degauss = {parameters.get('degauss', 0.01)}\n"
        "/\n"
    )

    # &ELECTRONS block
    electrons_block = (
        "&ELECTRONS\n"
        f"    mixing_beta = {parameters.get('mixing_beta', 0.7)}\n"
        "    conv_thr = 1.0e-10\n"
        "/\n"
    )

    # ATOMIC_SPECIES card
    atomic_species_lines = ["ATOMIC_SPECIES"]
    for symbol in unique_symbols:
        mass = atomic_masses[atoms.get_atomic_numbers()[atoms.get_chemical_symbols().index(symbol)]]
        pseudo = pseudopotentials[symbol]
        atomic_species_lines.append(f"{symbol}  {mass:.3f}  {pseudo}")
    atomic_species_card = "\n".join(atomic_species_lines) + "\n"

    # ATOMIC_POSITIONS card
    atomic_positions_lines = ["ATOMIC_POSITIONS {angstrom}"]
    for atom in atoms:
        atomic_positions_lines.append(
            f"{atom.symbol}   {atom.x:.8f}   {atom.y:.8f}   {atom.z:.8f}"
        )
    atomic_positions_card = "\n".join(atomic_positions_lines) + "\n"

    # K_POINTS card
    k_points_card = (
        "K_POINTS {automatic}\n"
        f"{kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0\n"
    )

    # CELL_PARAMETERS card
    cell_parameters_lines = ["CELL_PARAMETERS {angstrom}"]
    for vector in atoms.cell:
        cell_parameters_lines.append(
            f"  {vector[0]:.8f}   {vector[1]:.8f}   {vector[2]:.8f}"
        )
    cell_parameters_card = "\n".join(cell_parameters_lines) + "\n"

    return (
        f"{control_block}\n"
        f"{system_block}\n"
        f"{electrons_block}\n"
        f"{atomic_species_card}\n"
        f"{atomic_positions_card}\n"
        f"{k_points_card}\n"
        f"{cell_parameters_card}"
    )

def parse_qe_output(output: str) -> dict[str, Any]:
    """
    Parses the output of a Quantum Espresso calculation to extract key results.
    """
    if "convergence NOT achieved" in output:
        return {
            "was_successful": False,
            "error_message": "SCF did not converge",
        }

    try:
        total_energy_ry_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
        total_energy_ry = float(total_energy_ry_match.group(1))

        forces_block_match = re.search(
            r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n([\s\S]+?)\n\n", output
        )
        forces_str = forces_block_match.group(1)
        forces = []
        for line in forces_str.strip().split('\n'):
            parts = line.split()
            forces.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])

        stress_block_match = re.search(
            r"total stress.*?P=.*?\n(.*?)\n(.*?)\n(.*?)\n", output, re.DOTALL
        )
        stress = []
        for i in range(1, 4):
            line = stress_block_match.group(i).strip()
            parts = line.split()
            stress.append([float(parts[0]), float(parts[1]), float(parts[2])])

        return {
            "was_successful": True,
            "error_message": None,
            "total_energy_ev": total_energy_ry * RY_TO_EV,
            "forces": [[f * RY_TO_EV / 0.529177 for f in force] for force in forces],
            "stress": [[s * RY_BOHR3_TO_EV_ANG3 for s in row] for row in stress],
        }

    except (AttributeError, IndexError, ValueError) as e:
        return {
            "was_successful": False,
            "error_message": f"Parsing failed with exception: {e}",
        }
