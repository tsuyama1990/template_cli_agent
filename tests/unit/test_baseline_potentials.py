import numpy as np
from ase import Atoms
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential, calculate_zbl_potential

def test_calculate_lj_potential():
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [3.8, 0, 0]])
    epsilon = 0.0103  # eV
    sigma = 3.4  # Angstrom

    energy, forces = calculate_lj_potential(atoms, epsilon, sigma)

    # Check energy
    r = 3.8
    expected_energy = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    assert np.isclose(energy, expected_energy)

    # Check forces magnitude
    expected_force_mag = -4 * epsilon * (-12 * (sigma/r)**12 / r + 6 * (sigma/r)**6 / r)
    assert np.isclose(np.linalg.norm(forces[0]), abs(expected_force_mag))

def test_calculate_zbl_potential():
    atoms = Atoms('He2', positions=[[0, 0, 0], [1.0, 0, 0]])

    energy, forces = calculate_zbl_potential(atoms)

    assert isinstance(energy, float)
    assert energy > 0  # Repulsive potential
    assert isinstance(forces, np.ndarray)
    assert np.all(forces[0] == -forces[1])
