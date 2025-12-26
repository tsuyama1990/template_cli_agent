import pytest
from ase.build import bulk
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

def test_generate_qe_input():
    """
    Tests the generation of a Quantum Espresso input file.
    """
    atoms = bulk("Si", "diamond", a=5.43)
    qe_input = generate_qe_input(atoms)

    # Check for essential sections
    assert "&CONTROL" in qe_input
    assert "calculation = 'scf'" in qe_input
    assert "&SYSTEM" in qe_input
    assert "ibrav = 0" in qe_input
    assert "nat = 2" in qe_input
    assert "ntyp = 1" in qe_input
    assert "ecutwfc = 60.0" in qe_input
    assert "ATOMIC_SPECIES" in qe_input
    assert "Si  28.0855  Si.upf" in qe_input # Placeholder pseudopotential
    assert "ATOMIC_POSITIONS {angstrom}" in qe_input
    assert "CELL_PARAMETERS {angstrom}" in qe_input
    assert "K_POINTS {automatic}" in qe_input
    assert "6 6 6 0 0 0" in qe_input # Example k-points

def test_parse_qe_output_success():
    """
    Tests parsing of a successful Quantum Espresso output.
    """
    mock_output = """
    !    total energy              =     -150.0 Ry

    Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000001   -0.00000002    0.00000003
     atom    2 type  1   force =     0.00000001    0.00000002   -0.00000003

    total stress  (Ry/bohr**3)  (kbar)     P=   -0.00
    -0.00000001   0.00000000   0.00000000
     0.00000000  -0.00000001   0.00000000
     0.00000000   0.00000000  -0.00000001
    """
    result = parse_qe_output(mock_output)

    assert result.was_successful is True
    assert result.total_energy_ev == pytest.approx(-150.0 * 13.6057)
    assert len(result.forces) == 2
    assert result.forces[0][0] == pytest.approx(-0.00000001 * 51.4221)

def test_parse_qe_output_failure():
    """
    Tests parsing of a failed Quantum Espresso output.
    """
    mock_output = """
    Error in routine davcio (1):
     S matrix not positive definite
    """
    result = parse_qe_output(mock_output)

    assert result.was_successful is False
    assert "S matrix not positive definite" in result.error_message
