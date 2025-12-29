import numpy as np
import pytest
from ase.build import bulk

from mlip_autopipec.utils import dft_utils


@pytest.fixture
def si_atoms():
    """Fixture for a standard silicon Atoms object."""
    return bulk("Si", "diamond", a=5.43)

@pytest.fixture
def qe_parameters():
    """Fixture for a standard set of QE parameters."""
    return {
        "pseudopotentials": {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"},
        "k_points": [4, 4, 4],
        "ecutwfc": 60,
        "ecutrho": 240,
        "pseudo_dir": "./pseudos"
    }

def test_generate_qe_input_creates_valid_file(si_atoms, qe_parameters):
    """
    Tests that the generated QE input string contains all necessary sections and parameters.
    """
    input_str = dft_utils.generate_qe_input(si_atoms, qe_parameters)

    # Check for presence of key sections
    assert "&CONTROL" in input_str
    assert "&SYSTEM" in input_str
    assert "&ELECTRONS" in input_str
    assert "ATOMIC_SPECIES" in input_str
    assert "ATOMIC_POSITIONS" in input_str
    assert "K_POINTS" in input_str
    assert "CELL_PARAMETERS" in input_str

    # Check for specific parameter values
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 60" in input_str
    assert "ecutrho = 240" in input_str
    assert "Si.pbe-n-rrkjus_psl.1.0.0.UPF" in input_str
    assert "4 4 4 0 0 0" in input_str

def test_generate_qe_input_missing_params_raises_error(si_atoms):
    """
    Tests that a ValueError is raised if essential parameters are missing.
    """
    incomplete_params = {"k_points": [1, 1, 1]}
    with pytest.raises(ValueError, match="Missing one or more required parameters"):
        dft_utils.generate_qe_input(si_atoms, incomplete_params)

# Test data for parsing
SUCCESSFUL_QE_OUTPUT = """
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000000    0.00000000   -0.00000000
     atom    2 type  1   force =     0.00000000   -0.00000000    0.00000000

     !    total energy              =     -12.34567890 Ry

     total stress  (Ry/bohr**3)  (kbar)     P=      -0.12
       0.100000   0.000000   0.000000
       0.000000   0.100000   0.000000
       0.000000   0.000000   0.100000

     JOB DONE.
"""

FAILED_SCF_QE_OUTPUT = """
     convergence NOT achieved in  100 iterations

     End of self-consistent calculation
"""

NO_ENERGY_QE_OUTPUT = """
     Forces acting on atoms ...
     JOB DONE.
"""

def test_parse_qe_output_successful():
    """
    Tests parsing of a standard, successful QE output file.
    """
    result = dft_utils.parse_qe_output(SUCCESSFUL_QE_OUTPUT)

    assert result.was_successful is True
    assert result.error_message is None

    # Check energy ( -12.34567890 Ry * 13.60569... )
    assert pytest.approx(result.total_energy_ev, 0.001) == -167.971

    # Check forces ( should be near zero, converted to eV/A )
    assert len(result.forces) == 2
    assert np.allclose(result.forces, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    # Check stress ( 0.1 kbar * conversion factor )
    expected_stress_val = 0.1 / 160.21766208
    expected_stress = [[expected_stress_val, 0, 0], [0, expected_stress_val, 0], [0, 0, expected_stress_val]]
    assert len(result.stress) == 3
    assert np.allclose(result.stress, expected_stress)


def test_parse_qe_output_scf_failure():
    """
    Tests parsing of a QE output where the SCF cycle failed to converge.
    """
    result = dft_utils.parse_qe_output(FAILED_SCF_QE_OUTPUT)

    assert result.was_successful is False
    assert result.error_message == "SCF failed to converge."
    assert result.total_energy_ev == 0.0
    assert result.forces == []


def test_parse_qe_output_job_not_done():
    """
    Tests parsing when 'JOB DONE' is missing.
    """
    output = "This is an incomplete file."
    result = dft_utils.parse_qe_output(output)

    assert result.was_successful is False
    assert result.error_message == "Quantum Espresso job did not finish."


def test_parse_qe_output_missing_energy():
    """
    Tests that a parsing failure is correctly reported if energy is missing.
    """
    result = dft_utils.parse_qe_output(NO_ENERGY_QE_OUTPUT)

    assert result.was_successful is False
    assert "Failed to parse QE output: Could not find total energy." in result.error_message
