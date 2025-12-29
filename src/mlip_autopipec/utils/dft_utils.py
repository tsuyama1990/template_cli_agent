# ruff: noqa: D101, D102, D103, D104, D105, D107, C901
import re

import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.data.models import DFTResult


def generate_qe_input(
    atoms: Atoms,
    pseudo_potentials: dict[str, str],
    k_points: tuple[int, int, int] = (1, 1, 1),
    ecutwfc: float = 60.0,
) -> str:
    """Generates a Quantum Espresso input file string for a single-point calculation."""
    symbols = sorted(list(set(atoms.get_chemical_symbols())))
    from ase.data import atomic_masses, atomic_numbers

    species_block_lines = []
    for s in symbols:
        try:
            atomic_num = atomic_numbers[s]
            mass = atomic_masses[atomic_num]
            species_block_lines.append(f"  {s}  {mass:.4f}  {pseudo_potentials[s]}")
        except KeyError:
            raise ValueError(
                f"Atomic number or mass for symbol '{s}' not found in ase.data."
            )
    species_block = "\n".join(species_block_lines)

    positions_block = "\n".join(
        [f"  {atom.symbol}  {atom.x:.8f}  {atom.y:.8f}  {atom.z:.8f}" for atom in atoms]
    )

    cell_block = "\n".join(
        [f"    {row[0]:.8f}  {row[1]:.8f}  {row[2]:.8f}" for row in atoms.get_cell()]
    )

    input_str = f"""
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch',
    prefix = 'pwscf',
    outdir = './out',
    pseudo_dir = '.',
    verbosity = 'high',
    tprnfor = .true.,
    tstress = .true.,
/
&SYSTEM
    ibrav = 0,
    nat = {len(atoms)},
    ntyp = {len(symbols)},
    ecutwfc = {ecutwfc},
/
&ELECTRONS
    mixing_beta = 0.7,
    conv_thr = 1.0e-10,
/
ATOMIC_SPECIES
{species_block}
/
ATOMIC_POSITIONS {{angstrom}}
{positions_block}
/
K_POINTS {{automatic}}
  {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0
/
CELL_PARAMETERS {{angstrom}}
{cell_block}
/
"""
    return input_str


def parse_qe_output(output: str) -> DFTResult:
    """Parses the output of a Quantum Espresso calculation to extract results."""
    if "SCF NOT CONVERGED" in output:
        return DFTResult(
            was_successful=False,
            error_message="SCF failed to converge",
            total_energy_ev=0.0,
            forces=[],
            stress=[],
        )

    energy_match = re.search(r"!    total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
    if not energy_match:
        return DFTResult(
            was_successful=False,
            error_message="Could not find total energy in output.",
            total_energy_ev=0.0,
            forces=[],
            stress=[],
        )

    ry_to_ev = 13.605693122994
    total_energy = float(energy_match.group(1)) * ry_to_ev

    # Use re.findall for robustly capturing all force lines.
    force_lines = re.findall(r"force\s+=\s+(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)", output)

    # If no forces are found, check if the block header exists.
    # If not, the force calculation was likely not performed or failed.
    if not force_lines and "Forces acting on atoms" not in output:
        return DFTResult(
            was_successful=False,
            error_message="Could not find forces in output.",
            total_energy_ev=0.0,
            forces=[],
            stress=[],
        )

    ry_au_to_ev_a = 25.71104309541616 * 2
    forces = [list(map(float, line.split())) for line in force_lines]
    forces_ev_a = (np.array(forces) * ry_au_to_ev_a).tolist()

    stress_match = re.search(
        r"total stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=\s*(-?\d+\.\d+)\s*\n\s*"
        r"(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)",
        output,
    )

    if not stress_match:
        return DFTResult(
            was_successful=False,
            error_message="Could not find stress in output.",
            total_energy_ev=0.0,
            forces=[],
            stress=[],
        )

    ry_bohr3_to_ev_a3 = 14710.5073
    stress_values = list(map(float, stress_match.group(2).split()))
    stress_3x3 = np.array(
        [
            [stress_values[0], stress_values[5], stress_values[4]],
            [stress_values[5], stress_values[1], stress_values[3]],
            [stress_values[4], stress_values[3], stress_values[2]],
        ]
    )
    stress_ev_a3 = (stress_3x3 * ry_bohr3_to_ev_a3).tolist()

    return DFTResult(
        was_successful=True,
        total_energy_ev=total_energy,
        forces=forces_ev_a,
        stress=stress_ev_a3,
    )
