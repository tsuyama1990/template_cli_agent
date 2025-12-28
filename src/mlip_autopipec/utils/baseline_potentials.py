import numpy as np
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones


def get_lj_potential(atoms: Atoms) -> tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential for a given Atoms object.

    Args:
        atoms: The ASE Atoms object.

    Returns:
        A tuple containing the potential energy (float) and forces (np.ndarray).
    """
    calculator = LennardJones()
    atoms.set_calculator(calculator)

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    return energy, forces


def get_zbl_potential(atoms: Atoms) -> tuple[float, np.ndarray]:
    """
    Calculates the Ziegler-Biersack-Littmark (ZBL) universal potential.
    This is a simplified placeholder implementation for Cycle 01. A real
    implementation would be more complex.
    """
    # For Cycle 01, we will use the LJ potential as the baseline.
    return get_lj_potential(atoms)
