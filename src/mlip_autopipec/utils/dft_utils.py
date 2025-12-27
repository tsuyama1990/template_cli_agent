import re
import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_masses, atomic_numbers
from typing import Dict, Tuple

from ..data.models import DFTResult

# Conversion factors
RY_TO_EV = 13.605693009
RY_AU_TO_EV_A = RY_TO_EV / 0.52917721092
KBAR_TO_EV_A3 = 1.0 / 160.21766208

def generate_qe_input(
    atoms: Atoms,
    pseudo_dir: str,
    pseudopotentials: Dict[str, str],
    kpoints: Tuple[int, int, int],
    ecutwfc: float,
) -> str:
    """
    Generates a Quantum Espresso input file string for a single-point SCF calculation.
    """
    symbols = sorted(list(set(atoms.get_chemical_symbols())))

    # Basic input template
    input_str = f"""
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'pwscf'
    outdir = './out'
    pseudo_dir = '{pseudo_dir}'
    verbosity = 'high'
/
&SYSTEM
    ibrav = 0
    nat = {len(atoms)}
    ntyp = {len(symbols)}
    ecutwfc = {ecutwfc}
/
&ELECTRONS
    mixing_beta = 0.7
    conv_thr = 1.0e-8
/
ATOMIC_SPECIES
"""
    # Add atomic species
    for symbol in symbols:
        atomic_number = atomic_numbers[symbol]
        mass = atomic_masses[atomic_number]
        input_str += f"  {symbol}  {mass:.4f}  {pseudopotentials[symbol]}\n"

    # Add atomic positions
    input_str += "\nATOMIC_POSITIONS {angstrom}\n"
    for atom in atoms:
        input_str += f"  {atom.symbol} {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}\n"

    # Add cell parameters
    input_str += "\nCELL_PARAMETERS {angstrom}\n"
    for vector in atoms.get_cell():
        input_str += f"  {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}\n"

    # Add K-points
    input_str += "\nK_POINTS {automatic}\n"
    input_str += f"  {kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0\n"

    return input_str


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the output of a Quantum Espresso calculation to extract energy, forces, and stress.
    """
    if "not converged in" in output or "Maximum number of SCF cycles reached" in output:
        return DFTResult(
            total_energy_ev=0.0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message="SCF did not converge",
        )

    try:
        # Extract total energy
        energy_match = re.search(r"!\s+total energy\s+=\s+([-.\d]+)\s+Ry", output)
        if not energy_match:
            raise ValueError("Total energy not found in QE output.")
        total_energy_ry = float(energy_match.group(1))
        total_energy_ev = total_energy_ry * RY_TO_EV

        # Extract forces
        forces_section_match = re.search(
            r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n(.*?)\n\s*Total Force",
            output,
            re.DOTALL,
        )
        if not forces_section_match:
            raise ValueError("Forces not found in QE output.")

        forces_str = forces_section_match.group(1)
        forces = []
        for line in forces_str.strip().split('\n'):
            parts = line.split()
            force_ry_au = [float(f) for f in parts[-3:]]
            force_ev_a = [f * RY_AU_TO_EV_A for f in force_ry_au]
            forces.append(force_ev_a)

        # Extract stress
        stress_section_match = re.search(
            r"total stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=.*\n((?:\s*[-.\d]+\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+\n){3})",
            output,
        )
        if not stress_section_match:
            raise ValueError("Stress tensor not found in QE output.")

        stress_lines = stress_section_match.group(1).strip().split('\n')
        stress_kbar = np.array([list(map(float, line.split()[3:6])) for line in stress_lines])
        stress_ev_a3 = stress_kbar * KBAR_TO_EV_A3

        return DFTResult(
            total_energy_ev=total_energy_ev,
            forces=forces,
            stress=stress_ev_a3.tolist(),
            was_successful=True,
        )

    except (ValueError, IndexError) as e:
        return DFTResult(
            total_energy_ev=0.0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message=f"Failed to parse QE output: {e}",
        )
