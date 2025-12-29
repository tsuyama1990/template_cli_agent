# ruff: noqa: D101, D102, D103, D104, D105, D107

import numpy as np
from ase.atoms import Atoms


def calculate_lj_potential(
    atoms: Atoms, epsilon: float = 0.01032, sigma: float = 3.405
) -> tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential energy and forces for an Atoms object.

    Args:
        atoms: The ASE Atoms object.
        epsilon: LJ well depth in eV.
        sigma: LJ finite distance where inter-particle potential is zero, in Angstrom.

    Returns:
        A tuple containing the total potential energy (float) and forces (np.ndarray).
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)

    energy = 0.0
    forces = np.zeros((len(atoms), 3))

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            r = distances[i, j]

            if r == 0:
                continue

            sr6 = (sigma / r) ** 6
            sr12 = sr6**2

            energy += 4 * epsilon * (sr12 - sr6)

            # Force calculation
            # The force magnitude F = -dU/dr. A negative value indicates an attractive force.
            force_magnitude = 48 * epsilon / r * (sr12 - 0.5 * sr6)
            direction = (positions[j] - positions[i]) / r

            forces[i] -= force_magnitude * direction
            forces[j] += force_magnitude * direction

    return energy, forces
