from ase.build import bulk
from ase.data import atomic_masses, atomic_numbers

from mlip_autopipec.utils.dft_utils import generate_qe_input


def test_generate_qe_input_for_silicon():
    """Verify that the QE input generator produces a valid input for a silicon crystal."""
    # 1. Create a simple Silicon structure
    atoms = bulk("Si", "diamond", a=5.43)

    # 2. Call the function to generate the input string with specific k-points
    qe_input = generate_qe_input(
        atoms, calculation="scf", ecutwfc=40.0, kpts=(4, 4, 4)
    )

    # 3. Verify the output contains the expected sections and values
    assert "calculation = 'scf'" in qe_input
    assert "ecutwfc = 40.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input
    si_mass = f"{atomic_masses[atomic_numbers['Si']]:.4f}"
    assert f"    Si  {si_mass}  Si.pbe.UPF" in qe_input
    assert "ATOMIC_POSITIONS {crystal}" in qe_input
    assert "CELL_PARAMETERS {angstrom}" in qe_input
    assert "K_POINTS {automatic}" in qe_input
    assert "4 4 4 0 0 0" in qe_input


def test_generate_qe_input_defaults():
    """Verify that the function uses correct defaults for optional arguments."""
    atoms = bulk("Si", "diamond", a=5.43)
    qe_input = generate_qe_input(atoms)

    # Check that the default k-point mesh (Gamma point) is used
    assert "1 1 1 0 0 0" in qe_input
    # Check that the default pseudo_dir is './'
    assert "pseudo_dir = './'" in qe_input
