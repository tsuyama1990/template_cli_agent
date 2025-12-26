import pytest
import numpy as np
from ase import Atoms

from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential

def test_calculate_lj_potential():
    """
    Tests the Lennard-Jones potential calculation for a simple 2-atom system.
    """
    # Create a two-atom system (Argon-like parameters)
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 3.8]]) # Near equilibrium
    epsilon = 0.0103  # eV
    sigma = 3.4  # Angstrom

    energy, forces = calculate_lj_potential(atoms, epsilon=epsilon, sigma=sigma)

    # 1. Test energy calculation
    # At r = 3.8, r/sigma = 3.8/3.4 = 1.1176
    # V(r) = 4 * epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
    r_ratio = sigma / 3.8
    expected_energy = 4 * epsilon * (r_ratio**12 - r_ratio**6)
    assert energy == pytest.approx(expected_energy, abs=1e-5)

    # 2. Test force calculation
    # F(r) = -dV/dr = 48 * epsilon / r * [ (sigma/r)^12 - 0.5 * (sigma/r)^6 ]
    # Force is along the z-axis
    expected_force_mag = 48 * epsilon / 3.8 * (r_ratio**12 - 0.5 * r_ratio**6)

    # Force on atom 0 should be in the -z direction, force on atom 1 in +z
    assert forces[0][2] == pytest.approx(-expected_force_mag, abs=1e-5)
    assert forces[1][2] == pytest.approx(expected_force_mag, abs=1e-5)
    assert np.allclose(forces[0, :2], 0) # No force in x or y

    # 3. Test that forces sum to zero (Newton's 3rd law)
    assert np.allclose(np.sum(forces, axis=0), 0)

def test_lj_potential_cutoff():
    """Tests that the cutoff is respected."""
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 11.0]]) # Beyond cutoff
    energy, forces = calculate_lj_potential(atoms)

    assert energy == 0.0
    assert np.allclose(forces, 0.0)
