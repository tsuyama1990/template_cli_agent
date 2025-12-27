"""
Functions for calculating simple baseline interatomic potentials.
"""
import numpy as np
from ase.atoms import Atoms
from scipy.spatial.distance import pdist

def lennard_jones(atoms: Atoms, epsilon: float = 0.01032, sigma: float = 3.405) -> tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential energy and forces for an ASE Atoms object.

    Args:
        atoms: The ASE Atoms object.
        epsilon: The depth of the potential well in eV.
        sigma: The finite distance at which the inter-particle potential is zero, in Angstrom.

    Returns:
        A tuple containing:
        - total_energy (float): The total LJ potential energy in eV.
        - forces (np.ndarray): The forces on each atom in eV/Angstrom.
    """
    positions = atoms.get_positions()
    distances = pdist(positions)

    # Avoid division by zero for atoms at the same position
    distances = np.where(distances == 0, 1e-10, distances)

    r6 = (sigma / distances) ** 6
    r12 = r6 ** 2

    total_energy = np.sum(4 * epsilon * (r12 - r6))

    # Calculate forces
    forces = np.zeros_like(positions)
    num_atoms = len(atoms)

    if num_atoms < 2:
        return total_energy, forces

    # This is a naive O(n^2) implementation, sufficient for a baseline.
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            vec = positions[j] - positions[i]
            dist = np.linalg.norm(vec)
            if dist == 0:
                continue

            r_inv = 1.0 / dist
            s_r_6 = (sigma * r_inv) ** 6
            s_r_12 = s_r_6 ** 2

            force_magnitude = 24 * epsilon * r_inv * (2 * s_r_12 - s_r_6)
            force_vec = force_magnitude * (vec * r_inv)

            forces[i] += force_vec
            forces[j] -= force_vec

    return total_energy, forces
