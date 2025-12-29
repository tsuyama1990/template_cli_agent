# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential


def test_calculate_lj_potential_for_ar2():
    """
    Tests the Lennard-Jones calculation for a two-atom Argon system (Ar2 dimer).
    Uses standard LJ parameters for Argon.
    """
    # Standard parameters for Argon
    epsilon_ar = 0.01032  # eV
    sigma_ar = 3.405  # Angstrom
    equilibrium_dist = 2 ** (1 / 6) * sigma_ar

    # Create an Ar2 dimer at a distance slightly off equilibrium
    dist = equilibrium_dist * 1.1
    atoms = Atoms("Ar2", positions=[[0, 0, 0], [dist, 0, 0]])

    energy, forces = calculate_lj_potential(atoms, epsilon=epsilon_ar, sigma=sigma_ar)

    # --- Analytical Calculation for Verification ---
    r = dist
    sr6 = (sigma_ar / r) ** 6
    sr12 = sr6**2

    expected_energy = 4 * epsilon_ar * (sr12 - sr6)

    # The force is attractive. `dU/dr` is positive, so `F = -dU/dr` is negative.
    # This means the force vector on atom j should point in the opposite
    # direction of the interatomic vector (i->j)
    expected_force_magnitude = (
        -48 * epsilon_ar / r * (sr12 - 0.5 * sr6)
    )  # This is dU/dr, which is > 0
    expected_forces = np.array(
        [
            [expected_force_magnitude, 0.0, 0.0],  # Atom 0 is pulled in +x direction
            [-expected_force_magnitude, 0.0, 0.0],  # Atom 1 is pulled in -x direction
        ]
    )

    # Assert that the calculated values match the analytical ones
    assert np.isclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)
