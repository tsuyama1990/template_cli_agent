import numpy as np
from ase.atoms import Atoms
from typing import Tuple


def calculate_lj_potential(
    atoms: Atoms, epsilon: float = 0.0103, sigma: float = 3.4
) -> Tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential energy and forces for a set of atoms.
    Note: This is a very simple example and not a physically accurate baseline for most systems.

    Args:
        atoms: The ASE Atoms object.
        epsilon: The depth of the potential well (in eV).
        sigma: The distance at which the potential is zero (in Angstrom).

    Returns:
        A tuple containing the total energy (float) and forces (np.ndarray).
    """
    positions = atoms.get_positions()
    n_atoms = len(atoms)
    energy = 0.0
    forces = np.zeros((n_atoms, 3))

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist_vec = positions[i] - positions[j]
            dist = np.linalg.norm(dist_vec)

            if dist == 0:
                continue

            sr6 = (sigma / dist) ** 6
            sr12 = sr6**2

            energy += 4 * epsilon * (sr12 - sr6)

            # Force calculation: F = -dU/dr * (r_vec / r)
            force_mag = -24 * epsilon / dist * (2 * sr12 - sr6)
            force_vec = force_mag * (dist_vec / dist)
            forces[i] += force_vec
            forces[j] -= force_vec

    return energy, forces
