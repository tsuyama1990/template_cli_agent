import re

import numpy as np
import pytest

# This is the function we will implement in Phase 3
# We are writing the test for it now (TDD)
from mlip_autopipec.utils.dft_utils import parse_qe_output

# A realistic, but minimal, Quantum Espresso output for a simple Si calculation
SAMPLE_QE_OUTPUT = """
     Program PWSCF v.7.0 starts on  1Jan2024 at 12:00:00

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
     "P. Giannozzi et al., J. Phys.: Condens. Matter 21 395502 (2009);
     "P. Giannozzi et al., J. Phys.: Condens. Matter 29 465901 (2017);
     URL http://www.quantum-espresso.org",
     in publications or presentations arising from this work.

     Parallel version (MPI), running on     4 processors

     Waiting for input...
     Reading input from standard input

     ... (omitted setup and SCF cycles) ...

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000001   -0.00000001    0.00000001
     atom    2 type  1   force =     0.00000001    0.00000001   -0.00000001

     Total Force =     0.00000003     Total SCF correction =     0.00000000

     ...

!    total energy              =      -1.17122046 Ry

     ... (omitted timing information) ...

     The total stress is. (Ry/bohr**3)
     -0.00000000   -0.00000000   -0.00000000
     -0.00000000   -0.00000000   -0.00000000
     -0.00000000   -0.00000000   -0.00000000

     JOB DONE.
"""


def test_parse_qe_output():
    """
    Test the parsing of a sample Quantum Espresso output file.
    """
    parsed_data = parse_qe_output(SAMPLE_QE_OUTPUT)

    # Conversion factors from ASE
    RY_TO_EV = 13.605693122994
    RY_AU_TO_EV_ANG = 13.605693122994 / 0.529177210903

    # Check energy (with tolerance for floating point)
    expected_energy = -1.17122046 * RY_TO_EV
    assert parsed_data["energy"] == pytest.approx(expected_energy, abs=1e-6)

    # Check forces
    expected_forces = np.array([
        [-0.00000001, -0.00000001, 0.00000001],
        [0.00000001, 0.00000001, -0.00000001]
    ]) * RY_AU_TO_EV_ANG
    assert np.allclose(parsed_data["forces"], expected_forces, atol=1e-6)

    # Check stress (in this case, it's zero)
    # QE stress is in kbar, but ASE wants eV/Ang^3. The parser should handle this.
    # For now, we'll assume the parser returns it in a raw format that's easy to check.
    # The actual implementation will need to handle the conversion.
    # Since the raw stress is 0, the converted stress should also be 0.
    assert np.allclose(parsed_data["stress"], np.zeros((3, 3)), atol=1e-6)
