import pytest
import numpy as np
from ase import Atoms
from mlip_autopipec.utils.baseline_potentials import lennard_jones_potential

def test_lennard_jones_potential():
    """
    Tests the Lennard-Jones potential calculation for a two-atom system.
    """
    # Create a two-atom system at a known distance
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 3.8]]) # Approx. equilibrium

    # Lennard-Jones parameters for Argon
    epsilon = 0.0103  # eV
    sigma = 3.4  # Angstrom

    energy, forces = lennard_jones_potential(atoms, epsilon=epsilon, sigma=sigma)

    # 1. Verify Energy
    # At r=sigma, energy should be 0. We are slightly off.
    # U(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
    r = 3.8
    expected_energy = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    assert energy == pytest.approx(expected_energy, abs=1e-4)

    # 2. Verify Forces
    # Force is -dU/dr. It should be equal and opposite on the two atoms.
    assert np.allclose(forces[0], -forces[1])

    # Magnitude of the force on one atom along z-axis
    # F(r) = -d/dr [4*eps*((s/r)^12 - (s/r)^6)] = -4*eps*[-12*s^12/r^13 - (-6*s^6/r^7)]
    expected_force_mag = -4 * epsilon * (-12 * sigma**12 / r**13 + 6 * sigma**6 / r**7)
    assert forces[1][2] == pytest.approx(expected_force_mag)

def test_lennard_jones_no_interaction_for_single_atom():
    """
    Tests that a single-atom system has zero LJ energy and forces.
    """
    atoms = Atoms('Ar', positions=[[0, 0, 0]])
    energy, forces = lennard_jones_potential(atoms)

    assert energy == 0.0
    assert np.allclose(forces, [[0.0, 0.0, 0.0]])
