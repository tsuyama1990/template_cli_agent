import pytest
from ase.build import bulk
from ase.data import atomic_masses, atomic_numbers
from mlip_autopipec.utils.dft_utils import generate_qe_input

def test_generate_qe_input_for_silicon():
    """
    Tests the generation of a Quantum Espresso input file for a silicon crystal.
    This test is now robust and constructs the expected string programmatically,
    respecting the specific formatting required by Quantum Espresso.
    """
    # 1. Create a simple Silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # 2. Call the function to be tested
    qe_input = generate_qe_input(atoms, calculation='scf', ecutwfc=40.0)

    # 3. Verify the output against a programmatically generated, correct string
    assert isinstance(qe_input, str)
    assert "calculation = 'scf'" in qe_input
    assert "ecutwfc = 40.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input

    # --- Robust Assertion for QE Formatting ---
    # Construct the exact, correctly formatted species line that the function
    # SHOULD produce, matching the required column alignment and float precision.
    si_atomic_number = atomic_numbers['Si']
    si_mass = atomic_masses[si_atomic_number]
    expected_species_line = f"  {'Si':<4} {si_mass:10.4f}  Si.pbe.UPF"
    assert expected_species_line in qe_input
    # --- End Robust Assertion ---

    assert "ATOMIC_POSITIONS {crystal}" in qe_input
    assert "CELL_PARAMETERS {angstrom}" in qe_input
    assert "nat = 2" in qe_input
    assert "ntyp = 1" in qe_input
