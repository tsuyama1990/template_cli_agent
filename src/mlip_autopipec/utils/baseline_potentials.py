
import numpy as np
from ase import Atoms


def calculate_lj(atoms: Atoms, epsilon: float, sigma: float) -> tuple[float, np.ndarray]:
    """
    Calculates the total Lennard-Jones potential energy and forces for a given
    ASE Atoms object.

    Args:
        atoms: The ASE Atoms object.
        epsilon: The depth of the potential well in eV.
        sigma: The distance at which the potential is zero in Angstroms.

    Returns:
        A tuple containing:
        - The total potential energy in eV.
        - A NumPy array of the forces on each atom in eV/Angstrom.
    """
    positions = atoms.get_positions()
    n_atoms = len(atoms)

    total_energy = 0.0
    forces = np.zeros((n_atoms, 3))

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist_vec = positions[j] - positions[i]
            r = np.linalg.norm(dist_vec)

            if r == 0:
                continue

            sr = sigma / r
            sr6 = sr**6
            sr12 = sr6**2

            # Energy
            pair_energy = 4 * epsilon * (sr12 - sr6)
            total_energy += pair_energy

            # Force magnitude
            # F = -dE/dr = -4 * epsilon * [-12*(sigma^12/r^13) + 6*(sigma^6/r^7)]
            # F = (24 * epsilon / r) * [2 * (sigma/r)^12 - (sigma/r)^6]
            force_mag = (24 * epsilon / r) * (2 * sr12 - sr6)

            # Force vector
            force_vec = force_mag * (dist_vec / r)

            # Apply forces to the pair of atoms
            forces[i] += force_vec
            forces[j] -= force_vec

    return total_energy, forces
