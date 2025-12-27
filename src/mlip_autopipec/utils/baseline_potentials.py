import numpy as np
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones

def calculate_lj_potential(atoms: Atoms, epsilon: float = 0.01032, sigma: float = 3.405) -> tuple[float, np.ndarray]:
    """
    Calculates the energy and forces for an ASE Atoms object using the
    Lennard-Jones potential.

    This serves as a simple, fast baseline potential for the Delta Learning
    strategy. The default parameters are for Argon but provide a reasonable
    repulsive/attractive interaction for general-purpose use.

    Args:
        atoms: The ase.Atoms object for which to calculate the potential.
        epsilon: The depth of the potential well (in eV).
        sigma: The finite distance at which the inter-particle potential is zero (in Angstroms).

    Returns:
        A tuple containing:
        - The total potential energy (float, in eV).
        - The forces on each atom (np.ndarray, in eV/Angstrom).
    """
    # ASE's LennardJones calculator is a convenient and tested implementation.
    calculator = LennardJones(epsilon=epsilon, sigma=sigma, rc=10.0)
    atoms.calc = calculator

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    return energy, forces
