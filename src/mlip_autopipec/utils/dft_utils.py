import re
import subprocess
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_masses, atomic_numbers

from mlip_autopipec.data.models import DFTResult

# Conversion factors
RY_TO_EV = 13.605693122994
RY_AU_TO_EV_A = RY_TO_EV / 0.529177210903
RY_BOHR3_TO_EV_A3 = RY_TO_EV / (0.529177210903**3)


def generate_qe_input(
    atoms: Atoms,
    calculation: str = "scf",
    ecutwfc: float = 30.0,
    kpts: Optional[Tuple[int, int, int]] = None,
    pseudo_dir: str = "./",
    pseudopotentials: Optional[dict[str, str]] = None,
) -> str:
    """
    Generates a Quantum Espresso input string for the given atomic structure.

    Args:
        atoms: The ASE Atoms object.
        calculation: The type of calculation (e.g., 'scf', 'relax').
        ecutwfc: The wavefunction cutoff energy in Ry.
        kpts: The k-point mesh dimensions. Defaults to Gamma point (1, 1, 1).
        pseudo_dir: Path to the directory containing pseudopotential files.
        pseudopotentials: A dictionary mapping atomic symbols to pseudopotential filenames.
                          If None, it assumes a convention like 'Si' -> 'Si.pbe.UPF'.

    Returns:
        A string containing the formatted QE input file.
    """
    if kpts is None:
        kpts = (1, 1, 1)

    symbols = sorted(list(set(atoms.get_chemical_symbols())))
    if pseudopotentials is None:
        pseudopotentials = {s: f"{s}.pbe.UPF" for s in symbols}

    lines = []
    lines.append("&CONTROL")
    lines.append(f"    calculation = '{calculation}'")
    lines.append(f"    pseudo_dir = '{pseudo_dir}'")
    lines.append("    outdir = './out'")
    lines.append("/")
    lines.append("")

    lines.append("&SYSTEM")
    lines.append("    ibrav = 0")
    lines.append(f"    nat = {len(atoms)}")
    lines.append(f"    ntyp = {len(symbols)}")
    lines.append(f"    ecutwfc = {ecutwfc}")
    lines.append("/")
    lines.append("")

    lines.append("&ELECTRONS")
    lines.append("    conv_thr = 1.0e-10")
    lines.append("/")
    lines.append("")

    lines.append("ATOMIC_SPECIES")
    for symbol in symbols:
        mass = atomic_masses[atomic_numbers[symbol]]
        pseudo_file = pseudopotentials[symbol]
        lines.append(f"    {symbol}  {mass:.4f}  {pseudo_file}")
    lines.append("")

    lines.append("ATOMIC_POSITIONS {crystal}")
    scaled_positions = atoms.get_scaled_positions()
    for i, atom in enumerate(atoms):
        pos = scaled_positions[i]
        lines.append(f"    {atom.symbol}  {pos[0]:.9f}  {pos[1]:.9f}  {pos[2]:.9f}")
    lines.append("")

    lines.append("CELL_PARAMETERS {angstrom}")
    for vector in atoms.get_cell():
        lines.append(f"    {vector[0]:.9f}  {vector[1]:.9f}  {vector[2]:.9f}")
    lines.append("")

    lines.append("K_POINTS {automatic}")
    lines.append(f"    {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0")
    lines.append("")

    return "\n".join(lines)


def parse_qe_output(output: str) -> DFTResult:
    """
    Parses the standard output of a Quantum Espresso run to extract results.

    Args:
        output: The stdout string captured from the `pw.x` execution.

    Returns:
        A DFTResult object containing the parsed data. If the run failed,
        `was_successful` is False and `error_message` is populated.
    """
    if "convergence NOT achieved" in output or "ERROR" in output:
        return DFTResult(
            total_energy_ev=0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message=output[-500:],  # Capture the tail of the output for context
        )

    try:
        # Extract total energy
        energy_match = re.search(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", output)
        total_energy_ry = float(energy_match.group(1))

        # Extract forces
        forces_section = re.search(
            r"Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n([\s\S]+?)(?=total stress)",
            output,
        )
        forces_lines = forces_section.group(1).strip().split("\n")
        forces = []
        for line in forces_lines:
            parts = line.split()
            forces.append([float(parts[3]), float(parts[5]), float(parts[7])])

        # Extract stress
        stress_section = re.search(r"total stress\s+\(Ry/bohr\*\*3\).*?\n([\s\S]+)", output)
        stress_lines = stress_section.group(1).strip().split("\n")
        stress = []
        for line in stress_lines:
            # Filter out empty lines that might result from splitting
            if not line.strip():
                continue
            stress.append([float(p) for p in line.split()[:3]])

        return DFTResult(
            total_energy_ev=total_energy_ry * RY_TO_EV,
            forces=(np.array(forces) * RY_AU_TO_EV_A).tolist(),
            stress=(np.array(stress) * RY_BOHR3_TO_EV_A3).tolist(),
            was_successful=True,
        )
    except (AttributeError, IndexError, ValueError) as e:
        return DFTResult(
            total_energy_ev=0,
            forces=[],
            stress=[],
            was_successful=False,
            error_message=f"Parsing failed with error: {e}\nOutput tail:\n{output[-1000:]}",
        )


def run_qe(qe_command: str, input_file: str | Path) -> tuple[bool, str, str]:
    """
    Executes a Quantum Espresso calculation as a subprocess.

    Args:
        qe_command: The command to run QE (e.g., "pw.x", "mpirun -np 4 pw.x").
        input_file: The path to the QE input file.

    Returns:
        A tuple containing: (success_flag, stdout, stderr).
    """
    command = f"{qe_command} -in {input_file}"
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True, cwd=Path.cwd()
    )
    success = result.returncode == 0
    return success, result.stdout, result.stderr
