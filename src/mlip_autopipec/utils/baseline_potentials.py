import numpy as np
from ase.atoms import Atoms
from scipy.spatial.distance import pdist, squareform

def calculate_lj_potential(atoms: Atoms, epsilon: float = 0.0103, sigma: float = 3.4) -> tuple[float, np.ndarray]:
    """
    Calculates Lennard-Jones potential energy and forces.

    Returns:
        - Total energy (float) in eV.
        - Forces (np.ndarray) in eV/Angstrom.
    """
    positions = atoms.get_positions()
    distances = squareform(pdist(positions))

    # Avoid division by zero for self-distances
    distances[distances == 0] = np.inf

    # Efficient calculation of pairwise terms
    sigma_over_r = sigma / distances
    sigma_over_r6 = sigma_over_r**6
    sigma_over_r12 = sigma_over_r6**2

    # Energy calculation
    total_energy = np.sum(4.0 * epsilon * (sigma_over_r12 - sigma_over_r6)) / 2.0

    # Force calculation
    # F = 24 * epsilon * (2 * (sigma/r)^12 - (sigma/r)^6) / r^2 * r_vec
    force_magnitude = 24.0 * epsilon * (2.0 * sigma_over_r12 - sigma_over_r6) / (distances**2)

    forces = np.zeros_like(positions)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            r_vec = positions[i] - positions[j]
            force_vector = force_magnitude[i, j] * r_vec
            forces[i] += force_vector
            forces[j] -= force_vector

    return total_energy, forces

def calculate_zbl_potential(atoms: Atoms) -> tuple[float, np.ndarray]:
    """
    Calculates Ziegler-Biersack-Littmark (ZBL) universal repulsive potential.
    This is a simplified implementation. A full implementation would require
    element-specific parameters.

    Returns:
        - Total energy (float) in eV.
        - Forces (np.ndarray) in eV/Angstrom.
    """
    # Constants for a generic ZBL-like potential
    a = 0.46850  # Screening length
    e0 = 14.3996  # eV * Angstrom / e^2

    positions = atoms.get_positions()
    atomic_numbers = atoms.get_atomic_numbers()

    total_energy = 0.0
    forces = np.zeros_like(positions)

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            z1 = atomic_numbers[i]
            z2 = atomic_numbers[j]

            r_vec = positions[i] - positions[j]
            r = np.linalg.norm(r_vec)

            if r == 0:
                continue

            x = r / a

            # Screening function phi(x)
            phi = (0.1818 * np.exp(-3.2 * x) +
                   0.5099 * np.exp(-0.9423 * x) +
                   0.2802 * np.exp(-0.4028 * x) +
                   0.0281 * np.exp(-0.2016 * x))

            # Energy
            energy_ij = (e0 * z1 * z2 / r) * phi
            total_energy += energy_ij

            # Force derivative of phi(x)
            dphi_dx = (-3.2 * 0.1818 * np.exp(-3.2 * x) -
                       0.9423 * 0.5099 * np.exp(-0.9423 * x) -
                       0.4028 * 0.2802 * np.exp(-0.4028 * x) -
                       0.2016 * 0.0281 * np.exp(-0.2016 * x))

            dphi_dr = dphi_dx / a

            force_mag = -((e0 * z1 * z2 / r**2) * phi + (e0 * z1 * z2 / r) * dphi_dr)

            force_vec = force_mag * (r_vec / r)
            forces[i] += force_vec
            forces[j] -= force_vec

    return total_energy, forces
