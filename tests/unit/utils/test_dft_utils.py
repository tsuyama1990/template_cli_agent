from ase.build import bulk
import pytest
from pathlib import Path

from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

@pytest.fixture
def sample_qe_output():
    """Loads the sample QE output file."""
    path = Path(__file__).parent / "qe_output.txt"
    with open(path, "r") as f:
        return f.read()

def test_generate_qe_input():
    """Tests the generation of a Quantum Espresso input file."""
    atoms = bulk("Si", "diamond", a=5.43)
    pseudos = {"Si": "Si.upf"}
    kpts = (4, 4, 4)
    ecutwfc = 80.0

    input_str = generate_qe_input(atoms, pseudos, kpts, ecutwfc)

    assert "&SYSTEM" in input_str
    assert "ibrav = 0" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 80.0" in input_str
    assert "ATOMIC_SPECIES" in input_str
    assert "Si  1.0  Si.upf" in input_str
    assert "CELL_PARAMETERS angstrom" in input_str
    assert "ATOMIC_POSITIONS angstrom" in input_str
    assert "K_POINTS automatic" in input_str
    assert "4 4 4 0 0 0" in input_str
    assert "Si 0.000000000 0.000000000 0.000000000" in input_str

def test_parse_qe_output_successful(sample_qe_output):
    """Tests parsing a successful QE output."""
    result = parse_qe_output(sample_qe_output)

    assert result.was_successful is True
    assert result.error_message is None
    assert result.total_energy_ev == pytest.approx(-11.24355702 * 13.605693122994)
    assert len(result.forces) == 2
    assert len(result.forces[0]) == 3
    # QE output is in kbar, check one of the diagonal elements
    assert result.stress[0][0] == pytest.approx(-0.32)

def test_parse_qe_output_not_converged():
    """Tests parsing a QE output where SCF did not converge."""
    output = " JOB DONE." # missing convergence message
    result = parse_qe_output(output)
    assert result.was_successful is False
    assert "SCF did not converge" in result.error_message

def test_parse_qe_output_job_not_done():
    """Tests parsing an incomplete QE output."""
    output = "convergence has been achieved" # missing JOB DONE
    result = parse_qe_output(output)
    assert result.was_successful is False
    assert "did not finish" in result.error_message

def test_parse_qe_output_parsing_error():
    """Tests resilience against malformed output."""
    # This output has the convergence and JOB DONE flags, but malformed energy line
    output = "convergence has been achieved\n! total energy = garbage Ry\nForces acting on atoms ...\n JOB DONE."
    result = parse_qe_output(output)
    assert result.was_successful is False
    assert "Failed to parse output" in result.error_message
