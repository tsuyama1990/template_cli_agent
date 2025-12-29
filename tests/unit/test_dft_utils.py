import pytest
from ase.build import bulk
from mlip_autopipec.utils.dft_utils import generate_qe_input

def test_generate_qe_input_silicon_defaults():
    """
    Test case for generating a QE input for silicon with default settings.
    """
    silicon_atoms = bulk('Si', 'diamond', a=5.43)
    qe_input = generate_qe_input(silicon_atoms, calculation='scf', ecutwfc=40.0)

    # Verify standard sections
    assert "calculation = 'scf'" in qe_input
    assert "ecutwfc = 40.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input
    assert "Si " in qe_input
    assert "ATOMIC_POSITIONS" in qe_input
    assert "CELL_PARAMETERS" in qe_input

    # Verify K_POINTS section
    assert "K_POINTS {automatic}" in qe_input
    assert "4 4 4 0 0 0" in qe_input

    # Verify default pseudo_dir and outdir
    assert "pseudo_dir = './'" in qe_input
    assert "outdir = './'" in qe_input

    # Verify default pseudopotential
    assert "Si.pbe.UPF" in qe_input


def test_generate_qe_input_custom_args():
    """
    Test case for generating a QE input with custom arguments.
    """
    germanium_atoms = bulk('Ge', 'diamond', a=5.65)
    qe_input = generate_qe_input(
        germanium_atoms,
        k_points=(8, 8, 8),
        pseudo_dir='/tmp/pseudos',
        outdir='/tmp/output'
    )

    # Verify custom K_POINTS
    assert "K_POINTS {automatic}" in qe_input
    assert "8 8 8 0 0 0" in qe_input

    # Verify custom directories
    assert "pseudo_dir = '/tmp/pseudos'" in qe_input
    assert "outdir = '/tmp/output'" in qe_input


def test_generate_qe_input_custom_pseudos():
    """
    Test case for generating a QE input with custom pseudopotential filenames.
    """
    silicon_atoms = bulk('Si', 'diamond', a=5.43)
    custom_pseudos = {'Si': 'Si.custom.UPF'}
    qe_input = generate_qe_input(silicon_atoms, pseudos=custom_pseudos)

    # Verify custom pseudopotential is used
    assert "ATOMIC_SPECIES" in qe_input
    assert "Si.custom.UPF" in qe_input
    assert "Si.pbe.UPF" not in qe_input

def test_generate_qe_input_error_handling():
    """
    Test case for error handling in generate_qe_input.
    """
    # Pass a non-ASE object to trigger an error
    qe_input = generate_qe_input("not an atoms object")
    assert qe_input is None
