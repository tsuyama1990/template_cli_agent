"""Unit tests for the baseline_potentials module."""
import pytest
import numpy as np
from ase import Atoms

from mlip_autopipec.utils import baseline_potentials


def test_calculate_lennard_jones_argon_dimer():
    """
    Tests the Lennard-Jones calculation for a simple two-atom system (Argon).
    The equilibrium distance for LJ is 2^(1/6) * sigma.
    """
    # Parameters for Argon
    epsilon_ar = 0.0103  # eV
    sigma_ar = 3.4  # Angstrom
    equilibrium_dist = 2**(1/6) * sigma_ar

    # Create an argon dimer at the equilibrium distance
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [equilibrium_dist, 0, 0]])

    # Act
    total_energy, forces = baseline_potentials.calculate_lennard_jones(
        atoms, epsilon=epsilon_ar, sigma=sigma_ar
    )

    # Assert
    # At the equilibrium distance, the potential energy should be at its minimum, -epsilon.
    assert total_energy == pytest.approx(-epsilon_ar, abs=1e-6)

    # At the minimum of the potential, the forces should be zero.
    assert np.allclose(forces, np.zeros((2, 3)), atol=1e-6)

def test_calculate_lennard_jones_force_direction():
    """
    Tests that the force is repulsive when atoms are too close.
    """
    # Parameters for Argon
    epsilon_ar = 0.0103
    sigma_ar = 3.4
    compressed_dist = sigma_ar * 0.9  # Closer than the zero-potential distance

    # Create a compressed argon dimer
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [compressed_dist, 0, 0]])

    # Act
    _, forces = baseline_potentials.calculate_lennard_jones(
        atoms, epsilon=epsilon_ar, sigma=sigma_ar
    )

    # Assert
    # The force on the first atom should be in the negative x direction (repulsive).
    assert forces[0, 0] < 0
    # The force on the second atom should be in the positive x direction.
    assert forces[1, 0] > 0
    # The forces in y and z should be zero.
    assert np.allclose(forces[:, 1:], 0)
    # The total force on the system should be zero (Newton's third law).
    assert np.allclose(np.sum(forces, axis=0), 0)
