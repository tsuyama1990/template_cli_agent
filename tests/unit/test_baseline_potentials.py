import pytest
import numpy as np
from ase.atoms import Atoms
from src.mlip_autopipec.utils.baseline_potentials import lennard_jones_potential

def test_lennard_jones_potential():
    """
    Tests the Lennard-Jones potential calculation for a simple 2-atom system.
    """
    distance = 3.82
    # Create a simple, non-periodic two-atom system
    atoms = Atoms('Ar2', positions=[(0, 0, 0), (distance, 0, 0)])

    # Parameters for Argon
    epsilon = 0.0103  # eV
    sigma = 3.40  # Angstrom

    # Calculate potential
    energy, forces = lennard_jones_potential(atoms, epsilon=epsilon, sigma=sigma)

    # Analytical calculation for LJ potential energy and force magnitude
    r = distance
    s_over_r_6 = (sigma / r)**6
    expected_energy = 4 * epsilon * (s_over_r_6**2 - s_over_r_6)

    # Force is -dU/dr
    expected_force_magnitude = (24 * epsilon / r) * (2 * s_over_r_6**2 - s_over_r_6)

    # Total energy is per system, not per atom
    assert energy == pytest.approx(expected_energy)

    # Verify forces are equal and opposite
    assert np.allclose(forces[0], -forces[1])

    # Verify the magnitude of the force on one atom
    force_magnitude = np.linalg.norm(forces[0])
    assert force_magnitude == pytest.approx(abs(expected_force_magnitude))
