"""Utilities for interacting with Quantum Espresso."""
import re
from typing import Dict, Tuple, Optional, List

from ase.atoms import Atoms
from ase.data import atomic_masses, atomic_numbers
import numpy as np

from mlip_autopipec.data.models import DFTResult


def generate_qe_input(
    atoms: Atoms,
    pseudo_dir: str,
    pseudos: Dict[str, str],
    kpts: Tuple[int, int, int],
    ecutwfc: float,
    prefix: str = "qe",
    outdir: str = "out",
    calculation: str = "scf",
) -> str:
    """
    Generates a Quantum Espresso input file string from an ase.Atoms object.

    Args:
        atoms: The input atomic structure.
        pseudo_dir: Directory containing pseudopotential files.
        pseudos: A dictionary mapping atomic symbols to pseudopotential filenames.
        kpts: A 3-tuple specifying the k-point grid.
        ecutwfc: The plane-wave cutoff energy in Rydberg.
        prefix: A prefix for QE output files.
        outdir: The directory for QE output files.
        calculation: The type of QE calculation (e.g., 'scf', 'relax').

    Returns:
        A string containing the formatted Quantum Espresso input file.
    """
    symbols = sorted(list(set(atoms.get_chemical_symbols())))
    ntyp = len(symbols)
    nat = len(atoms)

    # Convert cell to Angstrom if it's not already
    cell_angstrom = atoms.get_cell()

    # CONTROL card
    control_card = f"""&CONTROL
    calculation = '{calculation}'
    restart_mode = 'from_scratch'
    prefix = '{prefix}'
    outdir = '{outdir}/'
    pseudo_dir = '{pseudo_dir}/'
    tstress = .true.
    tprnfor = .true.
/"""

    # SYSTEM card
    system_card = f"""&SYSTEM
    ibrav = 0
    nat = {nat}
    ntyp = {ntyp}
    ecutwfc = {ecutwfc}
/"""

    # ELECTRONS card
    electrons_card = """&ELECTRONS
    conv_thr = 1.0e-8
    mixing_beta = 0.7
/"""

    # CELL_PARAMETERS
    cell_card = "CELL_PARAMETERS {angstrom}\n"
    for vec in cell_angstrom:
        cell_card += f"  {vec[0]:.9f} {vec[1]:.9f} {vec[2]:.9f}\n"

    # ATOMIC_SPECIES
    species_card = "ATOMIC_SPECIES\n"
    for symbol in symbols:
        # Look up the correct atomic mass from ASE's database.
        atomic_num = atomic_numbers[symbol]
        mass = atomic_masses[atomic_num]
        species_card += f"  {symbol} {mass:.4f} {pseudos[symbol]}\n"

    # ATOMIC_POSITIONS
    positions_card = "ATOMIC_POSITIONS {angstrom}\n"
    for atom in atoms:
        pos = atom.position
        positions_card += f"  {atom.symbol} {pos[0]:.9f} {pos[1]:.9f} {pos[2]:.9f}\n"

    # K_POINTS
    kpoints_card = f"""K_POINTS {{automatic}}
  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0"""

    return "\n".join([
        control_card, system_card, electrons_card,
        cell_card, species_card, positions_card, kpoints_card
    ])


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the text output of a Quantum Espresso calculation.

    Args:
        output: The stdout from the 'pw.x' run.

    Returns:
        A DFTResult Pydantic model containing the parsed data.
    """
    # Check for specific, known failure modes first.
    if "convergence NOT achieved" in output:
        return DFTResult(was_successful=False, error_message="SCF calculation did not converge.")

    # Check for a clean exit.
    if "JOB DONE." not in output:
        return DFTResult(was_successful=False, error_message="Calculation did not finish (no 'JOB DONE.').")

    # Regex patterns
    energy_pattern = re.compile(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry")
    # This pattern is more robust: it captures from the header until a blank line, the energy block, or the stress block.
    forces_pattern = re.compile(r"Forces acting on atoms.*?\n(.*?)(?:\n\n|!    total energy|total stress)", re.DOTALL)
    stress_pattern = re.compile(r"total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*(-?\d+\.\d+)\n((?:\s+-?\d+\.\d+){6})\n")

    # Constants for unit conversion
    RY_TO_EV = 13.605693122994
    RY_BOHR_TO_EV_ANG = RY_TO_EV / 0.529177210903

    # Parse energy
    energy_match = energy_pattern.search(output)
    if not energy_match:
        return DFTResult(was_successful=False, error_message="Could not parse total energy.")
    total_energy_ev = float(energy_match.group(1)) * RY_TO_EV

    # Parse forces
    forces: Optional[List[List[float]]] = None
    forces_match = forces_pattern.search(output)
    if forces_match:
        forces_block = forces_match.group(1)
        force_lines = forces_block.strip().split('\n')
        parsed_forces = []
        for line in force_lines:
            parts = line.split()
            if len(parts) >= 9 and parts[0] == "atom":
                fx, fy, fz = map(float, parts[6:9])
                parsed_forces.append([fx * RY_BOHR_TO_EV_ANG, fy * RY_BOHR_TO_EV_ANG, fz * RY_BOHR_TO_EV_ANG])
        forces = parsed_forces

    # Parse stress
    stress: Optional[List[List[float]]] = None
    stress_match = stress_pattern.search(output)
    if stress_match:
        stress_voigt_kbar = np.fromstring(stress_match.group(2), sep=' ')
        # Convert kbar to eV/Ang^3
        stress_voigt_ev_ang3 = stress_voigt_kbar / 160.21766208
        # Convert Voigt (xx, yy, zz, xy, xz, yz) to 3x3 tensor
        s = stress_voigt_ev_ang3
        stress_tensor = [
            [s[0], s[5], s[4]],
            [s[5], s[1], s[3]],
            [s[4], s[3], s[2]]
        ]
        stress = stress_tensor

    return DFTResult(
        total_energy_ev=total_energy_ev,
        forces=forces,
        stress=stress,
        was_successful=True,
    )
