"""Unit tests for the DFT utilities."""
import pytest
from ase.atoms import Atoms
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output
import textwrap

@pytest.fixture
def sample_atoms():
    """Fixture for a sample Atoms object."""
    return Atoms('Si2',
                 cell=[[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]],
                 positions=[[0, 0, 0], [1.365, 1.365, 1.365]],
                 pbc=True)

@pytest.fixture
def sample_parameters():
    """Fixture for sample QE parameters."""
    return {
        "control": {
            "calculation": "'scf'",
            "prefix": "'silicon'",
            "pseudo_dir": "'./pseudos/'",
            "outdir": "'./out/'"
        },
        "system": {
            "ecutwfc": 60,
            "ecutrho": 240,
            "occupations": "'smearing'",
            "smearing": "'mv'",
            "degauss": 0.01
        },
        "electrons": {
            "conv_thr": "1.0e-10",
            "mixing_beta": 0.7
        },
        "k_points": {"scheme": "automatic", "grid": [4, 4, 4]}
    }

@pytest.fixture
def sample_pseudos():
    """Fixture for sample pseudopotentials."""
    return {"Si": "Si.upf"}


def test_generate_qe_input(sample_atoms, sample_parameters, sample_pseudos):
    """Test the generation of a QE input file."""
    input_str = generate_qe_input(sample_atoms, sample_parameters, sample_pseudos)

    assert "&CONTROL" in input_str
    assert "calculation = 'scf'" in input_str
    assert "&SYSTEM" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 60" in input_str
    assert "&ELECTRONS" in input_str
    assert "conv_thr = '1.0e-10'" in input_str
    assert "ATOMIC_SPECIES" in input_str
    assert "Si  1.0  Si.upf" in input_str
    assert "ATOMIC_POSITIONS {angstrom}" in input_str
    assert "Si  0.00000000  0.00000000  0.00000000" in input_str
    assert "Si  1.36500000  1.36500000  1.36500000" in input_str
    assert "K_POINTS {automatic}" in input_str
    assert "4 4 4 0 0 0" in input_str
    assert "CELL_PARAMETERS {angstrom}" in input_str
    assert "2.73000000  2.73000000  0.00000000" in input_str


def test_parse_qe_output_success():
    """Test parsing a successful QE output."""
    output = textwrap.dedent("""
        some preliminary text...
        !    total energy              =     -15.84501915 Ry
        some other text...
        Forces acting on atoms (cartesian axes, Ry/au):

         atom    1 type  1   force =    -0.00000133   -0.00000133   -0.00000133
         atom    2 type  1   force =     0.00000133    0.00000133    0.00000133

        Total force =     0.00000231

        total stress  (Ry/bohr**3)                (kbar)     P=       -1.54
        -0.00009852   0.00000000   0.00000000      -1.54   0.00   0.00
        0.00000000  -0.00009852   0.00000000       0.00  -1.54   0.00
        0.00000000   0.00000000  -0.00009852       0.00   0.00  -1.54

        JOB DONE.
    """)
    result = parse_qe_output(output)

    assert result.was_successful
    assert result.error_message is None
    assert pytest.approx(result.total_energy_ev, abs=1e-6) == -15.84501915 * 13.605693122994
    assert len(result.forces) == 2
    assert pytest.approx(result.forces[0][0], abs=1e-6) == -0.00000133 * 25.71104384
    assert len(result.stress) == 3
    ry_au3_to_ev_a3 = 13.605693122994 / (0.529177210903**3)
    assert pytest.approx(result.stress[0][0], abs=1e-6) == -0.00009852 * ry_au3_to_ev_a3

def test_parse_qe_output_scf_failure():
    """Test parsing a QE output where SCF failed."""
    output = textwrap.dedent("""
        ...
        convergence has not been achieved
        ...
    """)
    result = parse_qe_output(output)
    assert not result.was_successful
    assert result.error_message == "SCF did not converge"

def test_parse_qe_output_no_job_done():
    """Test parsing a QE output where JOB DONE is missing."""
    output = "some incomplete output"
    result = parse_qe_output(output)
    assert not result.was_successful
    assert "JOB DONE' not found" in result.error_message
