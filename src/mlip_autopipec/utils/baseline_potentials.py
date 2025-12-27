import numpy as np
from ase.atoms import Atoms
from typing import Tuple

def lennard_jones_potential(atoms: Atoms, epsilon: float, sigma: float) -> Tuple[float, np.ndarray]:
    """
    Calculates the total energy and forces for a system using the Lennard-Jones potential.

    Args:
        atoms: The ase.Atoms object.
        epsilon: The depth of the potential well (in eV).
        sigma: The distance at which the potential is zero (in Angstrom).

    Returns:
        A tuple containing:
        - The total potential energy of the system (float).
        - A numpy array of the forces on each atom.
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)
    forces = np.zeros_like(positions)
    total_energy = 0.0

    # Loop over unique pairs of atoms
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            r_vec = positions[j] - positions[i]

            # Minimum image convention for periodic systems
            cell = atoms.get_cell()
            if np.any(cell):
                r_vec -= cell @ np.round(np.linalg.inv(cell) @ r_vec)

            r = np.linalg.norm(r_vec)

            # Avoid division by zero
            if r == 0.0:
                continue

            s_over_r_6 = (sigma / r)**6

            # Energy calculation
            pair_energy = 4 * epsilon * (s_over_r_6**2 - s_over_r_6)
            total_energy += pair_energy

            # Force calculation (F = -dU/dr)
            force_magnitude = (24 * epsilon / r) * (2 * s_over_r_6**2 - s_over_r_6)
            force_vec = force_magnitude * (r_vec / r)

            forces[i] -= force_vec
            forces[j] += force_vec

    return total_energy, forces
