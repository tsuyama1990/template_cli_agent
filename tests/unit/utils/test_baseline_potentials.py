import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.utils import baseline_potentials


def test_calculate_lj_potential_argon_dimer():
    """
    Tests the LJ potential calculation for a simple Argon dimer.
    The LJ parameters for Argon are well-known (sigma=3.405 A, epsilon=0.0103 eV).
    ASE's default LJ calculator uses sigma=1.0, epsilon=1.0, which is not
    physically accurate but is sufficient to test that the function returns
    a deterministic, non-zero result.
    """
    # 1. Setup
    atoms = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 1.1]])

    # 2. Execute
    energy, forces = baseline_potentials.calculate_lj_potential(atoms)

    # 3. Verify
    # The expected energy for default ASE LJ params (sigma=1, eps=1) and r=1.1
    # is deterministic. The previous hand-calculation was incorrect.
    # The correct value is approximately -0.97789.
    assert isinstance(energy, float)
    assert pytest.approx(energy, 1e-5) == -0.97789300

    # Forces should be equal and opposite, pointing inwards (attractive)
    assert isinstance(forces, np.ndarray)
    assert forces.shape == (2, 3)
    assert np.allclose(forces[0], -forces[1]) # Equal and opposite
    assert forces[0][2] < 0  # Force on first atom is in -z direction (attractive)
    assert forces[1][2] > 0  # Force on second atom is in +z direction
