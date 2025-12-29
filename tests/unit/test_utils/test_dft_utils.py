# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
import re

import numpy as np
import pytest
from ase.atoms import Atoms

from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

# Fixtures with sample Quantum Espresso output strings


@pytest.fixture
def successful_qe_output() -> str:
    """Provides a sample of a successful QE pw.x output."""
    return """
     End of self-consistent calculation

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000011    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000011    0.00000000    0.00000000

     Total force =     0.00000016     Total SCF correction =     0.00000000

!    total energy              =     -78.82231388 Ry

     total stress  (Ry/bohr**3)                (kbar)     P=      -13.04
      0.00002138   -0.00000000   -0.00000000    -0.00000000   -0.00000000   -0.00000000

    """


@pytest.fixture
def failed_scf_qe_output() -> str:
    """Provides a sample of a QE output where SCF failed to converge."""
    return """
     iteration #  1     ecut=    60.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  4.0

     self-consistency has not been achieved: SCF NOT CONVERGED.
    """


@pytest.fixture
def atoms_h2() -> Atoms:
    """Provides a simple H2 Atoms object."""
    return Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)


def test_generate_qe_input(atoms_h2: Atoms):
    """Test the generation of a QE input file string."""
    pseudo_potentials = {"H": "H.pbe-rrkjus.UPF"}
    input_str = generate_qe_input(atoms_h2, pseudo_potentials)

    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "H  1.0080  H.pbe-rrkjus.UPF" in input_str
    assert "H  0.00000000  0.00000000  0.00000000" in input_str
    assert "H  0.00000000  0.00000000  0.74000000" in input_str
    assert "5.00000000  0.00000000  0.00000000" in input_str


def test_parse_qe_output_success(successful_qe_output: str):
    """Test parsing a successful QE output."""
    result = parse_qe_output(successful_qe_output)

    assert result.was_successful
    assert result.error_message is None
    assert np.isclose(result.total_energy_ev, -78.82231388 * 13.605693)
    assert len(result.forces) == 2
    assert len(result.forces[0]) == 3
    force_array = np.array([[-1.1e-7, 0, 0], [1.1e-7, 0, 0]]) * 25.711043 * 2
    assert np.allclose(np.array(result.forces), force_array, atol=1e-5)
    assert len(result.stress) == 3
    assert len(result.stress[0]) == 3


def test_parse_qe_output_scf_failure(failed_scf_qe_output: str):
    """Test parsing a QE output with an SCF convergence failure."""
    result = parse_qe_output(failed_scf_qe_output)

    assert not result.was_successful
    assert result.error_message == "SCF failed to converge"


def test_parse_qe_output_missing_energy():
    """Test parsing output where the energy line is missing."""
    output = "Some random text without energy."
    result = parse_qe_output(output)

    assert not result.was_successful
    assert result.error_message == "Could not find total energy in output."


def test_parse_qe_output_missing_forces(successful_qe_output: str):
    """Test parsing output where the forces block is missing."""
    # Simulate a case where forces are truly missing by removing the header and data
    output = successful_qe_output.replace("Forces acting on atoms", "Fcs acting on atoms")
    output = re.sub(r"force\s+=\s+.*", "", output)  # Corrected regex
    result = parse_qe_output(output)

    assert not result.was_successful
    assert result.error_message == "Could not find forces in output."


def test_parse_qe_output_missing_stress(successful_qe_output: str):
    """Test parsing output where the stress block is missing."""
    output = successful_qe_output.replace("total stress", "total pressure")
    result = parse_qe_output(output)

    assert not result.was_successful
    assert result.error_message == "Could not find stress in output."
