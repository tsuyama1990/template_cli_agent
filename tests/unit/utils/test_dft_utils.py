"""Unit tests for the dft_utils module."""
import pytest
from ase.build import bulk
from pathlib import Path

from mlip_autopipec.utils import dft_utils

# Sample data directory
SAMPLE_DATA_DIR = Path(__file__).parent / "sample_qe_outputs"
SAMPLE_DATA_DIR.mkdir(exist_ok=True)

# Create dummy sample output files for testing
SUCCESS_OUTPUT = """
Begin final coordinates
     atomic species an      at. pol. (x,y,z)             (a.u.)
 Si       1    ( 0.000000, 0.000000, 0.000000)
 Si       2    ( 2.715000, 2.715000, 2.715000)
End final coordinates

Forces acting on atoms (Ry/au):
     atom    1 type  1   force =    -0.00000000   -0.00000000   -0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000    0.00000000

!    total energy              =      -7.98393521 Ry

total stress  (Ry/bohr**3)                (kbar)     P=       -2.92
  -0.00001019   -0.00000000   -0.00000000    -0.00000000    0.00000000   -0.00000000

JOB DONE.
"""

CONVERGENCE_FAILURE_OUTPUT = """
     Maximum force       =     0.00000057    Threshold =     0.00000100
     Total force         =     0.00000057    Total SCF correction =     0.00000001
End of BFGS Geometry Optimization

convergence NOT achieved after  100 steps
JOB DONE.
"""

(SAMPLE_DATA_DIR / "success.out").write_text(SUCCESS_OUTPUT)
(SAMPLE_DATA_DIR / "failure.out").write_text(CONVERGENCE_FAILURE_OUTPUT)


def test_generate_qe_input():
    """Tests the generation of a QE input file for a simple Si cell."""
    si_atoms = bulk("Si", "diamond", a=5.43)
    pseudos = {"Si": "Si.upf"}
    input_str = dft_utils.generate_qe_input(
        atoms=si_atoms,
        pseudo_dir="/path/to/pseudos",
        pseudos=pseudos,
        kpts=(4, 4, 4),
        ecutwfc=60.0,
    )

    assert "&CONTROL" in input_str
    assert "calculation = 'scf'" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ecutwfc = 60.0" in input_str
    assert "ATOMIC_SPECIES" in input_str

    # Check mass with pytest.approx for float comparison robustness
    species_line = [line for line in input_str.split('\n') if line.strip().startswith("Si")][0]
    mass_str = species_line.strip().split()[1]
    # The atomic mass for Si in ASE's database is 28.085
    assert float(mass_str) == pytest.approx(28.085, abs=1e-4)
    assert "Si.upf" in species_line

    assert "ATOMIC_POSITIONS {angstrom}" in input_str
    assert "K_POINTS {automatic}" in input_str
    assert "4 4 4 0 0 0" in input_str


def test_parse_qe_output_success():
    """Tests parsing a successful QE output file."""
    output = (SAMPLE_DATA_DIR / "success.out").read_text()
    result = dft_utils.parse_qe_output(output)

    assert result.was_successful is True
    assert result.error_message is None
    assert result.total_energy_ev == pytest.approx(-7.98393521 * 13.605693, abs=1e-4)
    assert len(result.forces) == 2
    assert result.forces[0][0] == pytest.approx(0.0)
    assert len(result.stress) == 3
    assert result.stress[0][0] == pytest.approx(-0.00001019 / 160.21766, abs=1e-6)


def test_parse_qe_output_failure():
    """Tests parsing a failed QE output file."""
    output = (SAMPLE_DATA_DIR / "failure.out").read_text()
    result = dft_utils.parse_qe_output(output)

    assert result.was_successful is False
    assert result.error_message == "SCF calculation did not converge."

def test_parse_qe_no_job_done():
    """Tests parsing an incomplete QE output file."""
    output = "This is an incomplete file."
    result = dft_utils.parse_qe_output(output)

    assert result.was_successful is False
    assert result.error_message == "Calculation did not finish (no 'JOB DONE.')."
