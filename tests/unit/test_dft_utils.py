import pytest
from ase.build import bulk
from ase.data import atomic_masses, atomic_numbers

from mlip_autopipec.utils.dft_utils import generate_qe_input


def test_generate_qe_input_silicon():
    """
    Tests the generation of a Quantum Espresso input for a silicon crystal.
    """
    # 1. Create a simple Silicon structure
    atoms = bulk("Si", "diamond", a=5.43)

    # 2. Call generate_qe_input with specific k-points
    k_points = (4, 4, 4)
    qe_input = generate_qe_input(atoms, calculation="scf", ecutwfc=40.0, kpts=k_points)

    # 3. Verify the output
    assert "calculation = 'scf'" in qe_input
    assert "ecutwfc = 40.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input

    # Check for correct atomic mass for Si
    si_atomic_number = atomic_numbers["Si"]
    si_mass = atomic_masses[si_atomic_number]
    assert f"Si  {si_mass:.4f}  Si.pbe.UPF" in qe_input

    assert "ATOMIC_POSITIONS (crystal)" in qe_input

    # Check for correct k-points
    assert "K_POINTS automatic" in qe_input
    assert f"  {k_points[0]} {k_points[1]} {k_points[2]} 0 0 0" in qe_input
