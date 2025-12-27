import pytest
import numpy as np
from ase import Atoms
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential

def test_calculate_lj_potential_two_atoms():
    """Tests the LJ potential for a simple two-atom system."""
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 3.8]]) # A bit further than equilibrium
    epsilon = 0.0103  # eV for Argon
    sigma = 3.40  # Angstrom for Argon

    energy, forces = calculate_lj_potential(atoms, epsilon=epsilon, sigma=sigma)

    # Expected values from analytical calculation
    r = 3.8
    sr6 = (sigma / r)**6
    sr12 = sr6**2
    expected_energy = 4 * epsilon * (sr12 - sr6)
    expected_force_mag = -24 * epsilon * (2 * sr12 - sr6) / r

    assert np.isclose(energy, expected_energy)
    # Total force should be zero
    assert np.allclose(np.sum(forces, axis=0), 0)
    # Check magnitude of force on one atom
    assert np.isclose(np.linalg.norm(forces[0]), abs(expected_force_mag))

def test_calculate_lj_no_atoms():
    """Tests that the LJ potential returns zero for no atoms."""
    atoms = Atoms()
    energy, forces = calculate_lj_potential(atoms)
    assert energy == 0.0
    assert forces.shape == (0, 3)

def test_calculate_lj_one_atom():
    """Tests that the LJ potential returns zero for a single atom."""
    atoms = Atoms('Ar')
    energy, forces = calculate_lj_potential(atoms)
    assert energy == 0.0
    assert np.all(forces == 0.0)
