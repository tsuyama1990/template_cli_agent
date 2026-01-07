import re

import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTResult


def generate_qe_input(atoms: Atoms, parameters: dict) -> str:
    """
    Generates a Quantum Espresso input file string from an ASE Atoms object.
    """
    required = ['pseudopotentials', 'k_points', 'ecutwfc', 'ecutrho']
    if not all(p in parameters for p in required):
        raise ValueError(f"Missing one or more required parameters: {required}")

    lines = [
        "&CONTROL",
        "  calculation = 'scf'",
        "  restart_mode = 'from_scratch'",
        f"  pseudo_dir = '{parameters.get('pseudo_dir', '.')}/'",
        "  outdir = './out/'",
        "  prefix = 'qe'",
        "  tprnfor = .true.",
        "  tstress = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {len(atoms)}",
    ]

    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    lines.extend([
        f"  ntyp = {len(unique_symbols)}",
        f"  ecutwfc = {parameters['ecutwfc']}",
        f"  ecutrho = {parameters['ecutrho']}",
        "/",
        "&ELECTRONS",
        "  mixing_beta = 0.7",
        "  conv_thr = 1.0e-10",
        "/",
        "ATOMIC_SPECIES",
    ])

    for symbol in unique_symbols:
        pseudo = parameters['pseudopotentials'][symbol]
        mass = atoms[atoms.symbols.index(symbol)].mass
        lines.append(f"  {symbol} {mass:.4f} {pseudo}")

    lines.append("ATOMIC_POSITIONS {angstrom}")
    for atom in atoms:
        pos = atom.position
        lines.append(f"  {atom.symbol} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}")

    lines.append("K_POINTS {automatic}")
    k = parameters['k_points']
    lines.append(f"  {k[0]} {k[1]} {k[2]} 0 0 0")

    lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.get_cell():
        lines.append(f"  {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}")

    return "\n".join(lines)


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses Quantum Espresso output to extract energy, forces, and stress.
    """
    if "JOB DONE" not in output:
        error_msg = "Quantum Espresso job did not finish."
        if re.search(r'convergence NOT achieved', output):
            error_msg = "SCF failed to converge."
        return DFTResult(
            total_energy_ev=0.0, forces=[], stress=[], was_successful=False,
            error_message=error_msg
        )

    try:
        energy_match = re.search(r'!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry', output)
        if not energy_match:
            raise ValueError("Could not find total energy.")
        energy_ev = float(energy_match.group(1)) * 13.605693

        force_pattern = (
            r'Forces acting on atoms \(cartesian axes, Ry\/au\):\s*\n'
            r'(.*?)(?:\n\n|! total energy|total stress|\Z)'
        )
        forces_match = re.search(force_pattern, output, re.DOTALL)
        if not forces_match:
            raise ValueError("Could not find forces block.")

        force_lines = forces_match.group(1).strip().split('\n')
        forces_ry_au = [
            list(map(float, re.findall(r'[-]?\d+\.\d+', line)[-3:]))
            for line in force_lines if 'atom' in line
        ]
        forces_ev_a = np.array(forces_ry_au) * (13.605693 / 0.529177)

        stress_pattern = (
            r'total\s+stress\s+\(Ry\/bohr\*\*3\)\s*\(kbar\)\s*P=\s*[-]?\d+\.\d+\s*\n'
            r'\s*((?:[-]?\d+\.\d+\s*){3})\s*\n'
            r'\s*((?:[-]?\d+\.\d+\s*){3})\s*\n'
            r'\s*((?:[-]?\d+\.\d+\s*){3})'
        )
        stress_match = re.search(stress_pattern, output)
        if not stress_match:
            raise ValueError("Could not find stress block.")

        stress_lines = [stress_match.group(i).strip() for i in (1, 2, 3)]
        stress_kbar = np.array([list(map(float, line.split())) for line in stress_lines])
        stress_ev_a3 = stress_kbar / 160.21766

        return DFTResult(
            total_energy_ev=energy_ev,
            forces=forces_ev_a.tolist(),
            stress=stress_ev_a3.tolist(),
            was_successful=True,
        )

    except (ValueError, IndexError) as e:
        return DFTResult(
            total_energy_ev=0.0, forces=[], stress=[], was_successful=False,
            error_message=f"Failed to parse QE output: {e}"
        )
