"""Unit tests for the baseline potential utilities."""
import pytest
import numpy as np
from ase.atoms import Atoms
from mlip_autopipec.utils.baseline_potentials import lennard_jones

def test_lennard_jones_two_atoms():
    """Test the Lennard-Jones potential for a simple two-atom system."""
    # Place two atoms at a distance equal to sigma, where the potential should be zero
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [3.405, 0, 0]])
    energy, forces = lennard_jones(atoms)

    assert np.isclose(energy, 0.0)
    # At this distance, the force is attractive
    assert forces[0][0] > 0
    assert np.isclose(forces[0][0], -forces[1][0])

    # Place two atoms at the minimum of the potential, 2^(1/6) * sigma
    min_dist = 2**(1/6) * 3.405
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [min_dist, 0, 0]])
    energy, forces = lennard_jones(atoms)

    # Energy should be -epsilon
    assert np.isclose(energy, -0.01032)
    # Force should be zero at the minimum
    assert np.allclose(forces, 0.0)

def test_lennard_jones_single_atom():
    """Test that LJ potential returns zero for a single atom."""
    atoms = Atoms('Ar', positions=[[0, 0, 0]])
    energy, forces = lennard_jones(atoms)
    assert np.isclose(energy, 0.0)
    assert np.allclose(forces, 0.0)
