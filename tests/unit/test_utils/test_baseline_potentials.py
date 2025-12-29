import numpy as np
from ase import Atoms

from mlip_autopipec.utils import baseline_potentials


def test_lennard_jones_potential():
    """
    Tests the Lennard-Jones potential calculation for a simple dimer.
    Values are compared against a known analytical result.
    """
    # Parameters for Argon
    epsilon = 0.0103  # eV
    sigma = 3.4  # Angstrom

    distance = 3.8  # Angstrom, near the potential minimum
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [distance, 0, 0]])

    energy, forces = baseline_potentials.calculate_lj(
        atoms, epsilon=epsilon, sigma=sigma
    )

    # Analytical values
    r_sigma_6 = (sigma / distance) ** 6
    r_sigma_12 = r_sigma_6 ** 2
    expected_energy = 4 * epsilon * (r_sigma_12 - r_sigma_6)

    # Force is -dE/dr. The force points towards the other atom to reduce distance.
    expected_force_magnitude = -4 * epsilon * ((-12 * r_sigma_12 / distance) - (-6 * r_sigma_6 / distance))

    assert np.isclose(energy, expected_energy)

    # Check forces on both atoms
    assert np.isclose(forces[0][0], expected_force_magnitude)
    assert np.isclose(forces[1][0], -expected_force_magnitude)
    assert np.isclose(forces[0][1], 0.0) # No force in y or z
