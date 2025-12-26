from ase.atoms import Atoms
import numpy as np

def calculate_lj_potential(atoms: Atoms, epsilon: float = 0.0103, sigma: float = 3.4) -> tuple[float, np.ndarray]:
    """
    Calculates the Lennard-Jones potential energy and forces for a set of atoms.

    Args:
        atoms: The ASE Atoms object.
        epsilon: The depth of the potential well in eV.
        sigma: The distance at which the potential is zero in Angstroms.

    Returns:
        A tuple containing the total energy (float) and forces (np.ndarray).
    """
    positions = atoms.get_positions()
    distances = atoms.get_all_distances(mic=True)

    energy = 0.0
    forces = np.zeros((len(atoms), 3))

    # Get pairs of atoms to avoid double counting
    idx_i, idx_j = np.triu_indices(len(atoms), k=1)

    for i, j in zip(idx_i, idx_j):
        r = distances[i, j]
        if r > 10.0:  # Cutoff
            continue

        r6 = (sigma / r) ** 6
        r12 = r6 ** 2

        energy += 4 * epsilon * (r12 - r6)

        # Force calculation
        force_magnitude = 48 * epsilon / r * (r12 - 0.5 * r6)
        direction = (positions[i] - positions[j]) / r

        forces[i] += force_magnitude * direction
        forces[j] -= force_magnitude * direction

    return energy, forces
