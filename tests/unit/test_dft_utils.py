import pytest
from ase.build import bulk
from mlip_autopipec.utils.dft_utils import generate_qe_input

def test_generate_qe_input_for_silicon():
    """
    Tests the generation of a Quantum Espresso input file for a silicon crystal.
    Based on the UAT for Cycle 99 and feedback from code review.
    """
    # 1. Create a simple Silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # 2. Call generate_qe_input with explicit k-points
    qe_input = generate_qe_input(
        atoms,
        calculation='scf',
        ecutwfc=40.0,
        kpts=(6, 6, 6)
    )

    # 3. Verify the output
    assert "calculation = 'scf'" in qe_input
    assert "ecutwfc = 40.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input
    assert "Si  28.0850  Si.pbe.UPF" in qe_input
    assert "ATOMIC_POSITIONS {crystal}" in qe_input
    assert "CELL_PARAMETERS {angstrom}" in qe_input
    assert "ibrav = 0" in qe_input
    assert "nat = 2" in qe_input
    assert "ntyp = 1" in qe_input

    # Assert that the mandatory K_POINTS card is present and correct
    assert "K_POINTS {automatic}" in qe_input
    assert "6 6 6 0 0 0" in qe_input
