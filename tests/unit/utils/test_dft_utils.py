import pytest
import numpy as np
from ase import Atoms
import os
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

# Get the absolute path of the test directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TEST_DIR, '..', '..', 'data')

@pytest.fixture
def silicon_atoms():
    """Fixture for a simple silicon crystal."""
    return Atoms('Si2',
                 cell=[[0.0, 2.7, 2.7], [2.7, 0.0, 2.7], [2.7, 2.7, 0.0]],
                 positions=[[0, 0, 0], [1.35, 1.35, 1.35]],
                 pbc=True)

def test_generate_qe_input(silicon_atoms):
    """Tests the generation of a Quantum Espresso input file."""
    params = {
        'control': {'prefix': 'test'},
        'system': {'ecutwfc': 70.0},
    }
    pseudos = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

    input_str = generate_qe_input(silicon_atoms, params, pseudos)

    assert "&CONTROL" in input_str
    assert "prefix = 'test'" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 70.0" in input_str
    assert "ATOMIC_SPECIES" in input_str
    assert "Si " in input_str and "UPF" in input_str
    assert "CELL_PARAMETERS angstrom" in input_str
    assert "ATOMIC_POSITIONS angstrom" in input_str
    assert "K_POINTS automatic" in input_str

def test_parse_qe_output_success():
    """Tests parsing a successful Quantum Espresso output."""
    with open(os.path.join(DATA_DIR, 'qe_success.out'), 'r') as f:
        output = f.read()

    result = parse_qe_output(output)

    assert result.was_successful
    assert result.error_message is None
    assert np.isclose(result.total_energy_ev, -11.23032512 * 13.605693122994)
    assert len(result.forces) == 2
    assert np.allclose(result.forces, 0.0, atol=1e-8)
    assert np.isclose(result.stress[0][0], -0.00000001 * 13.605693122994 / (0.529177210903**3))

def test_parse_qe_output_failure():
    """Tests parsing a failed Quantum Espresso output."""
    with open(os.path.join(DATA_DIR, 'qe_fail.out'), 'r') as f:
        output = f.read()

    result = parse_qe_output(output)

    assert not result.was_successful
    assert "SCF did not converge" in result.error_message
    assert result.total_energy_ev == 0.0 # Should not parse an energy
    assert len(result.forces) == 0 # Should not parse forces

def test_parse_qe_output_job_done_but_no_energy():
    """Tests parsing an output where JOB DONE is present but energy is missing."""
    output = "This is a malformed output.\nJOB DONE."
    result = parse_qe_output(output)

    assert not result.was_successful
    assert "Could not parse energy or forces" in result.error_message
