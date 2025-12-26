import numpy as np
from ase.atoms import Atoms
from typing import Tuple

def lennard_jones_potential(
    atoms: Atoms, epsilon: float = 0.0103, sigma: float = 3.4
) -> Tuple[float, np.ndarray]:
    """
    Calculates the total energy and forces for a system using the Lennard-Jones potential.

    Args:
        atoms: The atomic structure.
        epsilon: The depth of the potential well in eV.
        sigma: The distance at which the potential is zero in Angstroms.

    Returns:
        A tuple containing:
        - The total potential energy in eV.
        - An array of forces on each atom in eV/Angstrom.
    """
    positions = atoms.get_positions()
    n_atoms = len(atoms)
    total_energy = 0.0
    forces = np.zeros((n_atoms, 3))

    if n_atoms < 2:
        return total_energy, forces

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Get the distance vector, respecting periodic boundary conditions
            dist_vec = atoms.get_distance(i, j, mic=True, vector=True)
            r = np.linalg.norm(dist_vec)

            if r < 1e-9:
                continue

            # Calculate LJ potential and force magnitude
            sr6 = (sigma / r) ** 6
            sr12 = sr6**2

            energy = 4 * epsilon * (sr12 - sr6)
            total_energy += energy

            # Force magnitude is -dU/dr
            force_mag = -4 * epsilon * (-12 * sr12 / r + 6 * sr6 / r)

            # Force vector
            unit_vec = dist_vec / r
            force_vec = force_mag * unit_vec

            # Apply forces according to Newton's third law
            forces[i] -= force_vec
            forces[j] += force_vec

    return total_energy, forces
