"""
Utilities for calculating simple, analytical baseline potentials.
"""
from typing import Tuple

import numpy as np
from ase.atoms import Atoms


def calculate_lennard_jones(
    atoms: Atoms, epsilon: float = 0.0103, sigma: float = 3.4
) -> Tuple[float, np.ndarray]:
    """
    Calculates the total energy and forces for a system using the Lennard-Jones potential.

    This is a simple pair potential. It is not intended to be accurate for real
    materials but serves as a baseline for the Delta Learning approach.

    Args:
        atoms: The ase.Atoms object for the system.
        epsilon: The depth of the potential well in eV.
        sigma: The distance at which the potential is zero in Angstroms.

    Returns:
        A tuple containing:
        - total_energy (float): The total potential energy of the system in eV.
        - forces (np.ndarray): A (N_atoms x 3) array of forces in eV/Angstrom.
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)

    # Exclude self-interaction
    np.fill_diagonal(distances, 0)

    # Avoid division by zero for zero distances
    r = np.where(distances == 0, 1e-10, distances)

    # Calculate LJ potential
    sigma_over_r6 = (sigma / r) ** 6
    sigma_over_r12 = sigma_over_r6 ** 2

    potential_energy_matrix = 4 * epsilon * (sigma_over_r12 - sigma_over_r6)
    np.fill_diagonal(potential_energy_matrix, 0)
    total_energy = np.sum(potential_energy_matrix) / 2.0  # Divide by 2 to avoid double counting

    # Calculate forces
    # Force = -dV/dr * (rij / |rij|)
    # dV/dr = 4 * epsilon * (-12 * sigma^12 / r^13 + 6 * sigma^6 / r^7)
    # The negative sign is applied here to get the force.
    force_magnitude = -4 * epsilon * (-12 * sigma_over_r12 / r + 6 * sigma_over_r6 / r)
    np.fill_diagonal(force_magnitude, 0)

    # Calculate force vectors
    forces = np.zeros_like(positions)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            r_vec = positions[j] - positions[i]
            # Handle periodic boundary conditions if necessary
            if np.any(atoms.pbc):
                dist, r_vec_pbc = atoms.get_distance(i, j, mic=True, vector=True)
                r_vec = r_vec_pbc

            # Normalize the direction vector
            r_norm = np.linalg.norm(r_vec)
            if r_norm > 0:
                direction = r_vec / r_norm
                force_ij = force_magnitude[i, j] * direction
                forces[i] -= force_ij
                forces[j] += force_ij

    return total_energy, forces
