import re
from io import StringIO
from typing import Any

import numpy as np
from ase import Atoms
from ase.io import write

from mlip_autopipec.data.models import DFTCompute

# Conversion factors
RY_TO_EV = 13.605693122994
RY_AU_TO_EV_A = RY_TO_EV / 0.529177210903
KBAR_TO_EV_A3 = 1 / 160.21766208


def create_qe_input_from_atoms(
    atoms: Atoms, config: DFTCompute, pseudopotentials: dict[str, str]
) -> str:
    """
    Creates a Quantum Espresso input file content from an ASE Atoms object.

    Args:
        atoms: The ASE Atoms object representing the structure.
        config: The DFTCompute configuration object.
        pseudopotentials: A dictionary mapping element symbols to their
                          pseudopotential filenames.

    Returns:
        A string containing the formatted Quantum Espresso input.
    """
    input_data = {
        'calculation': 'scf',
        'ecutwfc': config.ecutwfc,
        'ecutrho': config.ecutrho,
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.01,
        'tprnfor': True,
        'tstress': True,
    }

    density = int(round(config.kpoints_density))
    kpts = (density, density, density)

    string_io = StringIO()
    write(
        string_io,
        atoms,
        format='espresso-in',
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=kpts,
    )
    return string_io.getvalue()


def parse_qe_output(output_content: str) -> dict[str, Any] | None:
    """
    Parses the output of a Quantum Espresso (pw.x) calculation to extract
    energy, forces, and stress.

    Args:
        output_content: The full string content of the QE output file.

    Returns:
        A dictionary containing 'energy', 'forces', and 'stress' if successful,
        otherwise None.
    """
    try:
        energy = _parse_total_energy(output_content)
        forces = _parse_forces(output_content)
        stress = _parse_stress(output_content)

        if energy is None or forces is None or stress is None:
            return None

        return {'energy': energy, 'forces': forces, 'stress': stress}
    except (ValueError, IndexError):
        return None


def _parse_total_energy(content: str) -> float | None:
    """Parses the final total energy."""
    match = re.search(r'!\s+total energy\s+=\s+(-?[\d\.]+)\s+Ry', content)
    if match:
        return float(match.group(1)) * RY_TO_EV
    return None


def _parse_forces(content: str) -> np.ndarray | None:
    """Parses the forces on each atom."""
    force_block_match = re.search(
        r'Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n'
        r'(.*?)(?:\n\n|Total force|Writing forces)',
        content,
        re.DOTALL
    )
    if not force_block_match:
        return None

    force_block = force_block_match.group(1)
    lines = force_block.strip().split('\n')
    forces = []
    for line in lines:
        if 'atom' in line:
            parts = line.split('=')
            force_values = [float(f) for f in parts[-1].strip().split()]
            forces.append(force_values)

    if not forces:
        return None

    return np.array(forces) * RY_AU_TO_EV_A


def _parse_stress(content: str) -> np.ndarray | None:
    """Parses the total stress tensor and returns it in Voigt form."""
    stress_block_match = re.search(
        r'total stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*[-.\d]+\n'
        r'(.*?)(?=\n\n|\Z)',
        content,
        re.DOTALL
    )
    if not stress_block_match:
        return None

    stress_lines = stress_block_match.group(1).strip().split('\n')
    stress_matrix = np.array(
        [list(map(float, line.split()[:3])) for line in stress_lines]
    )

    voigt_stress = np.array([
        stress_matrix[0, 0],
        stress_matrix[1, 1],
        stress_matrix[2, 2],
        stress_matrix[1, 2],
        stress_matrix[0, 2],
        stress_matrix[0, 1]
    ])

    return voigt_stress * KBAR_TO_EV_A3
