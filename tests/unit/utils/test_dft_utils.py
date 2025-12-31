import pytest
from ase import Atoms
import numpy as np
import re

from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.utils.dft_utils import create_qe_input_from_atoms, parse_qe_output

# Fixture for a sample DFTCompute config
@pytest.fixture
def sample_dft_config():
    return DFTCompute(
        code="quantum_espresso",
        command="pw.x",
        pseudopotentials="SSSP_1.3_PBE_precision",
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=3.0
    )

# Fixture for a sample Atoms object
@pytest.fixture
def sample_atoms():
    return Atoms('Si', cell=[[2.7, 2.7, 0], [2.7, 0, 2.7], [0, 2.7, 2.7]], pbc=True, positions=[[0, 0, 0]])

# Sample QE output content is long, so it's defined once
SAMPLE_QE_OUTPUT = """
     End of self-consistent calculation

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000011    -0.00000011     0.00000022

     Total force =     0.00000027     Total SCF correction =     0.00000000

     !    total energy              =      -7.88151989 Ry

     total stress  (Ry/bohr**3)                (kbar)     P=      -1.12
      0.00007050   -0.00000000   -0.00000000   -0.00007050    0.00000000    0.00000000
     -0.00000000    0.00007050   -0.00000000    0.00000000   -0.00007050    0.00000000
     -0.00000000   -0.00000000    0.00007050    0.00000000    0.00000000   -0.00007050

"""

def test_create_qe_input_from_atoms_robust(sample_atoms, sample_dft_config):
    """
    Test the generation of a QE input file with robust parsing,
    ignoring whitespace and handling float precision.
    """
    pseudos = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}
    input_str = create_qe_input_from_atoms(sample_atoms, sample_dft_config, pseudos)

    def parse_namelist(content, name):
        match = re.search(rf'&{name.upper()}(.*?)/', content, re.DOTALL | re.IGNORECASE)
        assert match, f"Namelist '&{name}' not found."
        params = {}
        for line in match.group(1).strip().split('\n'):
            if '=' in line:
                parts = line.split('=')
                key = parts[0].strip()
                val = parts[1].strip().replace("'", "").replace(",", "")
                params[key] = val
        return params

    def parse_card(content, name):
        match = re.search(
            rf'^{name.upper()}.*?\n(.*?)(?=\n^[A-Z_]+|&|\Z)',
            content,
            re.DOTALL | re.IGNORECASE | re.MULTILINE
        )
        assert match, f"Card '{name}' not found."
        return [line.strip() for line in match.group(1).strip().split('\n') if line.strip()]

    # --- Verification ---
    control_params = parse_namelist(input_str, 'CONTROL')
    assert control_params['calculation'] == 'scf'

    system_params = parse_namelist(input_str, 'SYSTEM')
    assert float(system_params['ecutwfc']) == pytest.approx(60.0)
    assert float(system_params['ecutrho']) == pytest.approx(240.0)

    atomic_species = parse_card(input_str, 'ATOMIC_SPECIES')
    assert len(atomic_species) == 1
    species_parts = atomic_species[0].split()
    assert species_parts[0] == 'Si'
    assert float(species_parts[1]) == pytest.approx(28.085, abs=1e-2)
    assert species_parts[2] == 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'

    k_points = parse_card(input_str, 'K_POINTS')
    assert len(k_points) == 1
    assert " ".join(k_points[0].split()) == '3 3 3 0 0 0'

def test_parse_qe_output_success():
    """Test successful parsing of a QE output file."""
    results = parse_qe_output(SAMPLE_QE_OUTPUT)

    assert results is not None
    assert results['energy'] == pytest.approx(-7.88151989 * 13.605693122994)
    expected_forces = np.array([[-0.00000011, -0.00000011, 0.00000022]]) * (13.605693122994 / 0.529177210903)
    assert np.allclose(results['forces'], expected_forces)
    expected_stress_voigt = np.array([0.00007050, 0.00007050, 0.00007050, -0.0, -0.0, -0.0]) * (1 / 160.21766208)
    assert np.allclose(results['stress'], expected_stress_voigt)

def test_parse_qe_output_failure():
    """Test parsing of a failed/incomplete QE output."""
    results = parse_qe_output("This is not a valid QE output.")
    assert results is None

def test_parse_qe_output_no_energy():
    """Test parsing when energy is missing."""
    results = parse_qe_output(SAMPLE_QE_OUTPUT.replace("!    total energy", "some other text"))
    assert results is None

def test_parse_qe_output_no_forces():
    """Test parsing when forces are missing."""
    results = parse_qe_output(SAMPLE_QE_OUTPUT.replace("Forces acting on atoms", "Some other text"))
    assert results is None
