from ase import Atoms
from ase.data import atomic_masses
import numpy as np
from typing import Dict, Optional, Tuple

def _format_control_block(calculation: str, pseudo_dir: str, outdir: str) -> str:
    """Formats the &CONTROL block for the QE input file."""
    params = {
        'calculation': f"'{calculation}'",
        'pseudo_dir': f"'{pseudo_dir}'",
        'outdir': f"'{outdir}'",
    }
    lines = ["&CONTROL"]
    for key, value in params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")
    return "\n".join(lines) + "\n"

def _format_system_block(atoms: Atoms, ecutwfc: float) -> str:
    """Formats the &SYSTEM block for the QE input file."""
    if not atoms:
        raise ValueError("Input `atoms` object is empty.")
    params = {
        'ibrav': 0,
        'nat': len(atoms),
        'ntyp': len(set(atoms.get_chemical_symbols())),
        'ecutwfc': ecutwfc,
    }
    lines = ["&SYSTEM"]
    for key, value in params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")
    return "\n".join(lines) + "\n"

def _format_electrons_block() -> str:
    """Formats the &ELECTRONS block for the QE input file."""
    params = {'conv_thr': '1.0e-8'}
    lines = ["&ELECTRONS"]
    for key, value in params.items():
        lines.append(f"  {key} = {value}")
    lines.append("/")
    return "\n".join(lines) + "\n"

def _format_atomic_species(atoms: Atoms, pseudos: Optional[Dict[str, str]]) -> str:
    """Formats the ATOMIC_SPECIES card for the QE input file."""
    lines = ["ATOMIC_SPECIES"]
    unique_symbols = sorted(list(set(atoms.get_chemical_symbols())))
    atomic_numbers = atoms.get_atomic_numbers()
    symbols_list = np.array(atoms.get_chemical_symbols())
    for symbol in unique_symbols:
        atomic_number = atomic_numbers[np.where(symbols_list == symbol)[0][0]]
        mass = atomic_masses[atomic_number]
        pseudo_file = pseudos.get(symbol, f"{symbol}.pbe.UPF") if pseudos else f"{symbol}.pbe.UPF"
        lines.append(f"  {symbol} {mass:.4f} {pseudo_file}")
    return "\n".join(lines) + "\n"

def _format_atomic_positions(atoms: Atoms) -> str:
    """Formats the ATOMIC_POSITIONS card for the QE input file."""
    lines = ["ATOMIC_POSITIONS {crystal}"]
    scaled_positions = atoms.get_scaled_positions()
    symbols = atoms.get_chemical_symbols()
    for i in range(len(atoms)):
        lines.append(f"  {symbols[i]} {scaled_positions[i, 0]:.8f} {scaled_positions[i, 1]:.8f} {scaled_positions[i, 2]:.8f}")
    return "\n".join(lines) + "\n"

def _format_k_points(k_points: Tuple[int, int, int]) -> str:
    """Formats the K_POINTS card for the QE input file."""
    lines = ["K_POINTS {automatic}", f"  {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0"]
    return "\n".join(lines) + "\n"

def _format_cell_parameters(atoms: Atoms) -> str:
    """Formats the CELL_PARAMETERS card for the QE input file."""
    lines = ["CELL_PARAMETERS {angstrom}"]
    cell = atoms.get_cell()
    for i in range(3):
        lines.append(f"  {cell[i, 0]:.8f} {cell[i, 1]:.8f} {cell[i, 2]:.8f}")
    return "\n".join(lines)

def generate_qe_input(
    atoms: Atoms,
    calculation: str = 'scf',
    ecutwfc: float = 30.0,
    k_points: Tuple[int, int, int] = (4, 4, 4),
    pseudo_dir: str = './',
    outdir: str = './',
    pseudos: Optional[Dict[str, str]] = None
) -> Optional[str]:
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
        A string containing the formatted QE input file, or None if an error occurs.
    """
    try:
        input_parts = [
            _format_control_block(calculation, pseudo_dir, outdir),
            _format_system_block(atoms, ecutwfc),
            _format_electrons_block(),
            _format_atomic_species(atoms, pseudos),
            _format_atomic_positions(atoms),
            _format_k_points(k_points),
            _format_cell_parameters(atoms)
        ]
        return "\n".join(input_parts)
    except Exception as e:
        print(f"An error occurred in generate_qe_input: {e}")
        return None
