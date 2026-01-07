
import numpy as np
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones


def calculate_lj_potential(atoms: Atoms) -> tuple[float, np.ndarray]:
    """
    Calculates the potential energy and forces for a given Atoms object
    using the Lennard-Jones potential as implemented in ASE.

    Note: ASE's LJ calculator uses default parameters (sigma=1.0, epsilon=1.0)
    for element pairs if not specified. This is a simple baseline for the
    purpose of demonstrating Delta Learning. A real implementation might use
    more sophisticated or fitted parameters.

    Args:
        atoms: The ASE Atoms object.

    Returns:
        A tuple containing:
        - The total potential energy in eV.
        - The forces on each atom as a NumPy array.
    """
    # Detach any existing calculator to avoid interference
    atoms.calc = None

    # Use ASE's built-in LennardJones calculator
    calculator = LennardJones()
    atoms.calc = calculator

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    # Clean up by detaching the calculator
    atoms.calc = None

    return energy, forces
