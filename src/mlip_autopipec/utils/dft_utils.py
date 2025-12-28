from ase import Atoms
from ase.data import atomic_masses
import numpy as np
from typing import Dict, Optional, Tuple

def generate_qe_input(
    atoms: Atoms,
    calculation: str = 'scf',
    ecutwfc: float = 30.0,
    k_points: Tuple[int, int, int] = (4, 4, 4),
    pseudo_dir: str = './',
    outdir: str = './',
    pseudos: Optional[Dict[str, str]] = None
) -> str:
    """
    Generates a Quantum Espresso input string for the given atomic structure.

    Args:
        atoms: The ASE Atoms object.
        calculation: The type of calculation (e.g., 'scf', 'relax').
        ecutwfc: The wavefunction cutoff energy in Ry.
        k_points: A tuple of 3 integers for the k-point grid.
        pseudo_dir: The directory for pseudopotential files.
        outdir: The output directory.
        pseudos: A dictionary mapping atomic symbols to pseudopotential filenames.
                 If None, a default naming scheme is used (e.g., 'Si.pbe.UPF').

    Returns:
        A string containing the formatted QE input file.
    """

    control_params = {
        'calculation': f"'{calculation}'",
        'pseudo_dir': f"'{pseudo_dir}'",
        'outdir': f"'{outdir}'",
    }

    system_params = {
        'ibrav': 0,
        'nat': len(atoms),
        'ntyp': len(set(atoms.get_chemical_symbols())),
        'ecutwfc': ecutwfc,
    }

    electron_params = {
        'conv_thr': '1.0e-8',
    }

    # Create the input string
    input_str = ""

    # Control block
    input_str += "&CONTROL\n"
    for key, value in control_params.items():
        input_str += f"  {key} = {value}\n"
    input_str += "/\n\n"

    # System block
    input_str += "&SYSTEM\n"
    for key, value in system_params.items():
        input_str += f"  {key} = {value}\n"
    input_str += "/\n\n"

    # Electrons block
    input_str += "&ELECTRONS\n"
    for key, value in electron_params.items():
        input_str += f"  {key} = {value}\n"
    input_str += "/\n\n"

    # Atomic species
    input_str += "ATOMIC_SPECIES\n"
    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    atomic_numbers = atoms.get_atomic_numbers()
    symbols_list = np.array(atoms.get_chemical_symbols())
    for symbol in unique_symbols:
        atomic_number = atomic_numbers[np.where(symbols_list == symbol)[0][0]]
        mass = atomic_masses[atomic_number]
        if pseudos and symbol in pseudos:
            pseudo_file = pseudos[symbol]
        else:
            pseudo_file = f"{symbol}.pbe.UPF"
        input_str += f"  {symbol} {mass:.4f} {pseudo_file}\n"
    input_str += "\n"

    # Atomic positions
    input_str += "ATOMIC_POSITIONS {crystal}\n"
    scaled_positions = atoms.get_scaled_positions()
    symbols = atoms.get_chemical_symbols()
    for i in range(len(atoms)):
        input_str += f"  {symbols[i]} {scaled_positions[i, 0]:.8f} {scaled_positions[i, 1]:.8f} {scaled_positions[i, 2]:.8f}\n"
    input_str += "\n"

    # K-points
    input_str += "K_POINTS {automatic}\n"
    input_str += f"  {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0\n\n"

    # Cell parameters
    input_str += "CELL_PARAMETERS {angstrom}\n"
    cell = atoms.get_cell()
    for i in range(3):
        input_str += f"  {cell[i, 0]:.8f} {cell[i, 1]:.8f} {cell[i, 2]:.8f}\n"

    return input_str
