from ase.atoms import Atoms
from typing import Tuple
import re

from mlip_autopipec.data.models import DFTResult

def generate_qe_input(atoms: Atoms, pseudos: dict, kpts: Tuple[int, int, int], ecutwfc: float) -> str:
    """
    Generates a Quantum Espresso input file string for a given ASE Atoms object.
    """
    # This is a simplified input file generator. A real implementation would be
    # more robust, handling various options and pseudopotential formats.

    species = sorted(list(set(atoms.get_chemical_symbols())))

    input_lines = [
        "&CONTROL",
        "    calculation = 'scf'",
        "    restart_mode = 'from_scratch'",
        "    prefix = 'mlip'",
        "    outdir = './out'",
        "    pseudo_dir = './pseudos'",
        "    verbosity = 'high'",
        "/",
        "&SYSTEM",
        f"    ibrav = 0",
        f"    nat = {len(atoms)}",
        f"    ntyp = {len(species)}",
        f"    ecutwfc = {ecutwfc}",
        "/",
        "&ELECTRONS",
        "    mixing_beta = 0.7",
        "    conv_thr = 1.0e-8",
        "/",
        "ATOMIC_SPECIES",
    ]

    for s in species:
        input_lines.append(f"    {s}  1.0  {pseudos[s]}")

    input_lines.append("\nCELL_PARAMETERS angstrom")
    for vec in atoms.get_cell():
        input_lines.append(f"    {vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}")

    input_lines.append("\nATOMIC_POSITIONS angstrom")
    for atom in atoms:
        input_lines.append(f"    {atom.symbol} {atom.position[0]:.9f} {atom.position[1]:.9f} {atom.position[2]:.9f}")

    input_lines.append(f"\nK_POINTS automatic")
    input_lines.append(f"    {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0")

    return "\n".join(input_lines)


def parse_qe_output(stdout: str) -> DFTResult:
    """
    Parses the standard output of a Quantum Espresso calculation.
    """
    if "JOB DONE" not in stdout:
        return DFTResult(
            total_energy_ev=0.0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message="Calculation did not finish (no 'JOB DONE' found)."
        )

    # Check for convergence
    if "convergence has been achieved" not in stdout:
        return DFTResult(
            total_energy_ev=0.0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message="SCF did not converge."
        )

    try:
        # Extract total energy
        energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", stdout)
        energy_ry = float(energy_match.group(1))
        energy_ev = energy_ry * 13.605693122994

        # Extract forces
        forces_section_match = re.search(r"Forces acting on atoms \(Ry/au\):\s*\n(.*?)\n\s*Total force", stdout, re.DOTALL)
        forces_section = forces_section_match.group(1)
        force_lines = [line for line in forces_section.strip().split('\n') if 'atom' in line]
        forces = [list(map(float, line.split()[-3:])) for line in force_lines]
        # Convert forces from Ry/au to eV/Angstrom
        forces = [[f * 13.605693122994 / 0.529177210903 for f in force] for force in forces]


        # Extract stress tensor from the kbar column
        stress_section_match = re.search(r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*[\d\.\-]+\s*\n(.*?)\n", stdout, re.DOTALL)
        stress_section = stress_section_match.group(1)
        stress_lines = stress_section.strip().split('\n')
        stress_kbar = [list(map(float, line.split()[3:6])) for line in stress_lines]


        return DFTResult(
            total_energy_ev=energy_ev,
            forces=forces,
            stress=stress_kbar,
            was_successful=True
        )

    except (AttributeError, IndexError, ValueError) as e:
        return DFTResult(
            total_energy_ev=0.0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message=f"Failed to parse output: {e}"
        )
