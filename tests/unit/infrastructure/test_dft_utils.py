"""Unit tests for DFT utility functions."""

import pytest
import numpy as np
from mlip_autopipec.infrastructure.dft_utils import parse_qe_output, QEOutputError

# Sample QE output for a successful run (simplified)
SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =     -11.45567302 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000   -0.00000000

total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
 0.00000004   -0.00000000    0.00000000
 0.00000000    0.00000004    0.00000000
 0.00000000    0.00000000    0.00000004
"""

SAMPLE_FAILURE_OUTPUT = "convergence NOT achieved"

def test_parse_qe_output_success():
    """Test successful parsing of QE output."""
    result = parse_qe_output(SAMPLE_SUCCESS_OUTPUT)
    assert result.energy == -11.45567302

def test_parse_qe_output_failure():
    """Test that parsing a failed run raises an exception."""
    with pytest.raises(QEOutputError, match="SCF calculation did not converge."):
        parse_qe_output(SAMPLE_FAILURE_OUTPUT)
