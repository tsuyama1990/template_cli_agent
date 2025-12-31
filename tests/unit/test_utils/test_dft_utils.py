"""Unit tests for DFT utility functions."""

import numpy as np
import pytest

from mlip_autopipec.utils.dft_utils import QEOutputError, parse_qe_output

# Sample QE output for a successful run (simplified)
SAMPLE_SUCCESS_OUTPUT = """
     Program PWSCF v.6.5 starts on 12Jan2024 at 10:00:00

     ... some QE header ...

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000   -0.00000000
     atom    2   force =     0.00000101    0.00000000    0.00000000

     ... some other output ...

     !    total energy              =     -11.45567302 Ry

     ... more output ...

     total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
      0.00000004   -0.00000000    0.00000000   -0.00000001    0.00000000   -0.00000000
      0.00000000    0.00000004    0.00000000    0.00000000   -0.00000001    0.00000000
      0.00000000    0.00000000    0.00000004    0.00000000    0.00000000   -0.00000001

"""

# Sample QE output for a failed (non-converged) run
SAMPLE_FAILURE_OUTPUT = """
     Program PWSCF v.6.5 starts on 12Jan2024 at 10:05:00

     ... some QE header ...

     iteration #  1     ecut=    40.00 Ry     beta= 0.70
     total cpu time spent up to now is        0.1 secs

     End of self-consistent calculation

     convergence NOT achieved after  100 iterations: stopping

"""


def test_parse_qe_output_success():
    """Test successful parsing of QE output."""
    result = parse_qe_output(SAMPLE_SUCCESS_OUTPUT)

    assert isinstance(result, object)
    assert result.energy == -11.45567302

    expected_forces = np.array([[-0.00000101, 0.0, 0.0], [0.00000101, 0.0, 0.0]])
    np.testing.assert_allclose(result.forces, expected_forces, atol=1e-8)

    expected_stress = np.array(
        [[0.00000004, -0.0, 0.0], [0.0, 0.00000004, 0.0], [0.0, 0.0, 0.00000004]]
    )
    np.testing.assert_allclose(result.stress, expected_stress, atol=1e-8)


def test_parse_qe_output_failure():
    """Test that parsing a failed run raises an exception."""
    with pytest.raises(QEOutputError, match="SCF calculation did not converge."):
        parse_qe_output(SAMPLE_FAILURE_OUTPUT)


def test_parse_qe_output_missing_energy():
    """Test handling of output missing the energy block."""
    bad_output = SAMPLE_SUCCESS_OUTPUT.replace("!    total energy", "!    something else")
    with pytest.raises(QEOutputError, match="Could not find total energy"):
        parse_qe_output(bad_output)


def test_parse_qe_output_missing_forces():
    """Test handling of output missing the forces block."""
    bad_output = SAMPLE_SUCCESS_OUTPUT.replace("Forces acting on atoms", "Something else")
    with pytest.raises(QEOutputError, match="Could not find forces"):
        parse_qe_output(bad_output)


def test_parse_qe_output_missing_stress():
    """Test handling of output missing the stress block."""
    bad_output = SAMPLE_SUCCESS_OUTPUT.replace("total stress", "something else")
    with pytest.raises(QEOutputError, match="Could not find stress tensor"):
        parse_qe_output(bad_output)
