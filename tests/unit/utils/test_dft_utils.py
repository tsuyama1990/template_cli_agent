import re

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.utils.dft_utils import (
    KBAR_TO_EV_A3,
    RY_AU_TO_EV_A,
    RY_TO_EV,
    create_qe_input_from_atoms,
    parse_qe_output,
)


@pytest.fixture
def sample_dft_config():
    return DFTCompute(
        code="quantum_espresso",
        command="pw.x",
        pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'},
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=3.0
    )

@pytest.fixture
def sample_atoms():
    return Atoms(
        'Si',
        cell=[[2.7, 2.7, 0], [2.7, 0, 2.7], [0, 2.7, 2.7]],
        pbc=True,
        positions=[[0, 0, 0]]
    )

SAMPLE_QE_OUTPUT = """
!    total energy              =      -7.88151989 Ry
Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000011    -0.00000011     0.00000022

total stress  (Ry/bohr**3)                (kbar)     P=      -1.12
  0.00007050   -0.00000000   -0.00000000
 -0.00000000    0.00007050   -0.00000000
 -0.00000000   -0.00000000    0.00007050
"""

def test_create_qe_input_from_atoms_robust(sample_atoms, sample_dft_config):
    pseudos = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}
    input_str = create_qe_input_from_atoms(sample_atoms, sample_dft_config, pseudos)

    def parse_namelist(content, name):
        match = re.search(
            rf'&{name.upper()}(.*?)/', content, re.DOTALL | re.IGNORECASE
        )
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
        return [
            line.strip() for line in match.group(1).strip().split('\n') if line.strip()
        ]

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
    results = parse_qe_output(SAMPLE_QE_OUTPUT)

    assert results is not None
    assert results['energy'] == pytest.approx(-7.88151989 * RY_TO_EV)
    expected_forces = np.array(
        [[-0.00000011, -0.00000011, 0.00000022]]
    ) * RY_AU_TO_EV_A
    assert np.allclose(results['forces'], expected_forces)
    expected_stress = np.array(
        [0.00007050, 0.00007050, 0.00007050, 0.0, 0.0, 0.0]
    ) * KBAR_TO_EV_A3
    assert np.allclose(results['stress'], expected_stress)

def test_parse_qe_output_failure():
    results = parse_qe_output("This is not a valid QE output.")
    assert results is None

def test_parse_qe_output_no_energy():
    output = SAMPLE_QE_OUTPUT.replace("!    total energy", "some other text")
    results = parse_qe_output(output)
    assert results is None

def test_parse_qe_output_no_forces():
    output = SAMPLE_QE_OUTPUT.replace("Forces acting on atoms", "Some other text")
    results = parse_qe_output(output)
    assert results is None
