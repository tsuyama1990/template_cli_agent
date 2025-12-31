from ase import Atoms
from ase.io import write
from io import StringIO
import re
import numpy as np
from typing import Dict, Any, Optional, Tuple

from mlip_autopipec.data.models import DFTCompute

def create_qe_input_from_atoms(atoms: Atoms, config: DFTCompute, pseudopotentials: Dict[str, str]) -> str:
    """
    Creates a Quantum Espresso input file content from an ASE Atoms object.

    Args:
        atoms: The ASE Atoms object representing the structure.
        config: The DFTCompute configuration object.
        pseudopotentials: A dictionary mapping element symbols to their pseudopotential filenames.

    Returns:
        A string containing the formatted Quantum Espresso input.
    """
    # ASE's writer requires a dictionary of parameters.
    # We construct this from our Pydantic model.
    input_data = {
        'calculation': 'scf',  # For single point energy/force calculation
        'ecutwfc': config.ecutwfc,
        'ecutrho': config.ecutrho,
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.01,
        'tprnfor': True, # To print forces
        'tstress': True, # To print stress
    }

    # K-points density should be converted to an integer grid
    density = int(round(config.kpoints_density))
    kpts = (density, density, density)

    # Use StringIO to capture the output of ase.io.write as a string
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


def parse_qe_output(output_content: str) -> Optional[Dict[str, Any]]:
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

def _parse_total_energy(content: str) -> Optional[float]:
    """Parses the final total energy."""
    match = re.search(r'!\s+total energy\s+=\s+(-?[\d\.]+)\s+Ry', content)
    if match:
        # Convert from Rydberg to eV
        return float(match.group(1)) * 13.605693122994
    return None

def _parse_forces(content: str) -> Optional[np.ndarray]:
    """Parses the forces on each atom."""
    force_block_match = re.search(
        r'Forces acting on atoms \(cartesian axes, Ry/au\):\s*\n(.*?)(?:\n\n|Total force|Writing forces)',
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

    # Convert from Ry/au to eV/Angstrom
    return np.array(forces) * (13.605693122994 / 0.529177210903)

def _parse_stress(content: str) -> Optional[np.ndarray]:
    """Parses the total stress tensor and returns it in Voigt form."""
    stress_block_match = re.search(
        r'total stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*[-.\d]+\n(.*?)(?=\n\n|\Z)',
        content,
        re.DOTALL
    )
    if not stress_block_match:
        return None

    stress_lines = stress_block_match.group(1).strip().split('\n')
    stress_matrix = np.array([list(map(float, line.split()[:3])) for line in stress_lines])

    # ASE expects stress in Voigt order: [xx, yy, zz, yz, xz, xy]
    # QE output is in matrix form:
    # xx xy xz
    # yx yy yz
    # zx zy zz
    # The stress tensor is symmetric (xy=yx, etc.)
    voigt_stress = np.array([
        stress_matrix[0, 0],
        stress_matrix[1, 1],
        stress_matrix[2, 2],
        stress_matrix[1, 2],
        stress_matrix[0, 2],
        stress_matrix[0, 1]
    ])

    # Convert from kbar to eV/Angstrom^3
    return voigt_stress * (1 / 160.21766208)
