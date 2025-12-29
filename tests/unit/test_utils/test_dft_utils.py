import pytest
from ase import Atoms

from mlip_autopipec.utils import dft_utils


def test_generate_qe_input_creates_correct_string_for_simple_atom():
    """
    Tests that the generated QE input string is syntactically correct and
    contains the right values for a simple, single-atom system.
    """
    atoms = Atoms('H', positions=[[0, 0, 0]], cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]], pbc=True)

    # These would typically come from a config file
    pseudopotentials = {'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'}
    parameters = {
        'calculation': 'scf',
        'ecutwfc': 40,
        'ecutrho': 160,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,
        'mixing_beta': 0.7,
    }
    kpoints = (1, 1, 1)

    expected_control_block = (
        "&CONTROL\n"
        "    calculation = 'scf'\n"
        "    restart_mode = 'from_scratch'\n"
        "    prefix = 'pwscf'\n"
        "    outdir = './out'\n"
        "    pseudo_dir = './pseudos'\n"
        "/\n"
    )

    expected_system_block = (
        "&SYSTEM\n"
        "    ibrav = 0\n"
        "    nat = 1\n"
        "    ntyp = 1\n"
        "    ecutwfc = 40\n"
        "    ecutrho = 160\n"
        "    occupations = 'smearing'\n"
        "    smearing = 'gaussian'\n"
        "    degauss = 0.01\n"
        "/\n"
    )

    expected_electrons_block = (
        "&ELECTRONS\n"
        "    mixing_beta = 0.7\n"
        "    conv_thr = 1.0e-10\n"
        "/\n"
    )

    expected_atomic_species = (
        "ATOMIC_SPECIES\n"
        "H  1.008  H.pbe-rrkjus_psl.1.0.0.UPF\n"
    )

    expected_atomic_positions = (
        "ATOMIC_POSITIONS {angstrom}\n"
        "H   0.00000000   0.00000000   0.00000000\n"
    )

    expected_kpoints = (
        "K_POINTS {automatic}\n"
        "1 1 1 0 0 0\n"
    )

    expected_cell_parameters = (
        "CELL_PARAMETERS {angstrom}\n"
        "  10.00000000   0.00000000   0.00000000\n"
        "  0.00000000   10.00000000   0.00000000\n"
        "  0.00000000   0.00000000   10.00000000\n"
    )

    qe_input = dft_utils.generate_qe_input(atoms, parameters, pseudopotentials, kpoints)

    # Check for presence and basic correctness of each block
    assert expected_control_block in qe_input
    assert expected_system_block in qe_input
    assert expected_electrons_block in qe_input
    assert expected_atomic_species in qe_input
    assert expected_atomic_positions in qe_input
    assert expected_kpoints in qe_input
    assert expected_cell_parameters in qe_input


def test_parse_qe_output_success():
    """Tests parsing a successful QE run output."""
    with open("tests/unit/test_data/qe_outputs/h_atom_successful_run.out") as f:
        output = f.read()

    result = dft_utils.parse_qe_output(output)

    assert result["was_successful"] is True
    assert result["error_message"] is None
    assert pytest.approx(result["total_energy_ev"], 1e-6) == -0.45206260 * 13.6057
    assert len(result["forces"]) == 1
    assert len(result["forces"][0]) == 3
    assert pytest.approx(result["forces"][0][0], 1e-6) == 0.0
    # Check the stress conversion from Ry/bohr^3 to eV/A^3
    expected_stress_ev_ang3 = -0.00000005 * dft_utils.RY_BOHR3_TO_EV_ANG3
    assert pytest.approx(result["stress"][0][0], 1e-6) == expected_stress_ev_ang3


def test_parse_qe_output_failure_convergence():
    """Tests parsing a failed (non-converged) QE run output."""
    with open("tests/unit/test_data/qe_outputs/h_atom_failed_run.out") as f:
        output = f.read()

    result = dft_utils.parse_qe_output(output)

    assert result["was_successful"] is False
    assert "SCF did not converge" in result["error_message"]
    assert "total_energy_ev" not in result
