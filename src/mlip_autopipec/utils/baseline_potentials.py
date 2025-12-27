import numpy as np
from ase.atoms import Atoms
from typing import Tuple

def calculate_lj_potential(atoms: Atoms, epsilon: float = 0.01, sigma: float = 3.4) -> Tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential energy and forces for a set of atoms.

    Args:
        atoms: The ASE Atoms object.
        epsilon: The depth of the potential well in eV.
        sigma: The distance at which the potential is zero in Angstroms.

    Returns:
        A tuple containing the total energy (float) and the forces (np.ndarray).
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)

    energy = 0.0
    forces = np.zeros((len(atoms), 3))

    # Get pairs of atoms to avoid double counting
    idx_i, idx_j = np.triu_indices(len(atoms), k=1)

    for i, j in zip(idx_i, idx_j):
        dist = distances[i, j]
        if dist == 0: continue # Should not happen with triu_indices k=1

        # Calculate LJ potential
        sr6 = (sigma / dist) ** 6
        sr12 = sr6 ** 2
        energy += 4 * epsilon * (sr12 - sr6)

        # Calculate force magnitude
        force_mag = -24 * epsilon * (2 * sr12 - sr6) / dist

        # Calculate force vector
        direction = (positions[i] - positions[j]) / dist
        force_vec = force_mag * direction

        forces[i] += force_vec
        forces[j] -= force_vec

    return energy, forces
