import re
from ase.atoms import Atoms
from typing import Dict, Tuple, Optional

# SSSP pseudopotential recommendations for simplicity
SSSP_PSEUDOS = {
    'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
    'H': 'H_pbe_v1.4.uspp.F.UPF',
    'O': 'O.pbe-nl-rrkjus_psl.1.0.0.UPF',
}
SSSP_CUTOFFS = {
    'Si': (60, 240),
    'H': (50, 200),
    'O': (60, 240),
}

def generate_qe_input(atoms: Atoms, pseudos: Dict[str, str], ecutwfc: int, ecutrho: int, k_points: Tuple[int, int, int]) -> str:
    """
    Generates a Quantum Espresso input file string for a given ASE Atoms object.

    Args:
        atoms: The ASE Atoms object for the calculation.
        pseudos: A dictionary mapping element symbols to pseudopotential filenames.
        ecutwfc: The kinetic energy cutoff for wavefunctions (Rydbergs).
        ecutrho: The kinetic energy cutoff for charge density (Rydbergs).
        k_points: A tuple of 3 integers for the k-point mesh.

    Returns:
        A string containing the formatted Quantum Espresso input file.
    """
    atom_types = sorted(list(set(atoms.get_chemical_symbols())))
    nat = len(atoms)
    ntyp = len(atom_types)

    input_str = f"""
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'pwscf'
    outdir = './out'
    pseudo_dir = './pseudos'
/
&SYSTEM
    ibrav = 0
    nat = {nat}
    ntyp = {ntyp}
    ecutwfc = {ecutwfc}
    ecutrho = {ecutrho}
/
&ELECTRONS
    mixing_beta = 0.7
    conv_thr = 1.0e-8
/
ATOMIC_SPECIES
"""
    for symbol in atom_types:
        # A placeholder mass; QE doesn't use it for SCF calculations
        mass = 1.0
        pseudo_filename = pseudos[symbol]
        input_str += f"    {symbol}  {mass:.4f} {pseudo_filename}\n"

    input_str += "\nCELL_PARAMETERS angstrom\n"
    for vec in atoms.cell:
        input_str += f"    {vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}\n"

    input_str += "\nATOMIC_POSITIONS angstrom\n"
    for atom in atoms:
        input_str += f"    {atom.symbol} {atom.position[0]:.9f} {atom.position[1]:.9f} {atom.position[2]:.9f}\n"

    input_str += f"\nK_POINTS automatic\n"
    input_str += f"    {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0\n"

    return input_str.strip()

def parse_qe_output(output: str) -> Tuple[Optional[float], Optional[list], Optional[list], bool, Optional[str]]:
    """
    Parses the standard output of a Quantum Espresso calculation.

    Args:
        output: The string captured from stdout.

    Returns:
        A tuple containing:
        - Total energy in eV (or None if not found).
        - Forces as a list of lists (or None).
        - Stress tensor (or None).
        - A boolean indicating if the calculation was successful.
        - An error message string (or None).
    """
    energy = None
    forces = None
    stress = None
    error_message = None

    if "JOB DONE" not in output:
        # Look for a common convergence failure message
        if "SCF NOT CONVERGED" in output:
            error_message = "SCF did not converge."
        else:
            error_message = "QE job did not complete successfully. Check logs."
        return None, None, None, False, error_message

    # Find total energy
    energy_match = re.search(r'!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry', output)
    if energy_match:
        # Convert from Rydberg to eV
        energy = float(energy_match.group(1)) * 13.6057

    # Find forces
    force_match = re.search(r'Forces acting on atoms \(cartesian axes, Ry\/au\):\s*\n\s*\n((?:\s*atom\s+\d+\s+type\s+\w+\s+force\s+=\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n)+)', output)
    if force_match:
        forces = []
        force_lines = force_match.group(1).strip().split('\n')
        for line in force_lines:
            parts = line.split()
            # Convert from Ry/au to eV/Angstrom
            fx = float(parts[6]) * 51.422067
            fy = float(parts[7]) * 51.422067
            fz = float(parts[8]) * 51.422067
            forces.append([fx, fy, fz])

    # Find stress tensor
    stress_match = re.search(r'total\s+stress\s+\(Ry\/bohr\*\*3\)\s+pbc\s+=\s*\n((?:\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n)+)', output, re.MULTILINE)
    if stress_match:
        stress = []
        stress_lines = stress_match.group(1).strip().split('\n')
        # Convert from Ry/bohr^3 to eV/Angstrom^3
        conversion_factor = 13.6057 / (0.529177**3)
        for line in stress_lines:
            parts = line.split()
            s_xx = float(parts[0]) * conversion_factor
            s_yy = float(parts[1]) * conversion_factor
            s_zz = float(parts[2]) * conversion_factor
            stress.append([s_xx, s_yy, s_zz])
            # Note: QE outputs the tensor in a different order than ASE expects
            # This is a simplified parsing; a real implementation would be more robust.

    # Simple check for success
    was_successful = energy is not None and forces is not None

    return energy, forces, stress, was_successful, error_message
