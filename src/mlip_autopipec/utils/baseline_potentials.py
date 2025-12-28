# Description: Implementation of simple analytical baseline potentials for Delta Learning.
import numpy as np
from ase.atoms import Atoms


def zbl_potential(atoms: Atoms) -> tuple[float, np.ndarray]:
    """
    Calculates the energy and forces using the Ziegler-Biersack-Littmark (ZBL) potential.

    This is a simplified implementation for demonstration purposes. It uses a generic
    screening function and does not depend on the specific element types.

    Args:
        atoms: The ASE Atoms object for the structure.

    Returns:
        A tuple containing:
        - The total potential energy (float).
        - The forces on each atom (np.ndarray of shape (N, 3)).
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)
    natoms = len(atoms)

    energy = 0.0
    forces = np.zeros((natoms, 3))
    k_e = 14.399645  # Coulomb's constant in eV * Angstrom / e^2

    # A generic Z value for demonstration. A real implementation would use atomic numbers.
    z = 14.0

    for i in range(natoms):
        for j in range(i + 1, natoms):
            r_ij = distances[i, j]
            if r_ij > 1e-6:  # Avoid division by zero
                # Generic screening function phi(x) = 0.5 * exp(-0.3 * x)
                # where x = r / a_s and a_s is the screening length.
                # Simplified here for demonstration.
                a_s = 0.2
                x = r_ij / a_s
                phi = 0.5 * np.exp(-0.3 * x)

                # Energy
                pair_energy = (k_e * z**2 / r_ij) * phi
                energy += pair_energy

                # Force
                dphi_dr = -0.5 * (0.3 / a_s) * np.exp(-0.3 * x)
                force_magnitude = (
                    -k_e * z**2 * (phi / r_ij**2 - dphi_dr / r_ij)
                )

                # Vector pointing from j to i
                r_vec = positions[i] - positions[j]
                # Minimum image convention adjustment
                cell = atoms.get_cell()
                r_vec -= cell @ np.round(np.linalg.solve(cell.T, r_vec.T)).T

                force_vec = force_magnitude * (r_vec / r_ij)
                forces[i] += force_vec
                forces[j] -= force_vec

    return energy, forces
