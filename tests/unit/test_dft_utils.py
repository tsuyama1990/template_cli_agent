# tests/unit/test_dft_utils.py

import numpy as np
import pytest

from mlip_autopipec.utils.dft_utils import parse_qe_output

# Sample QE output for a 2-atom Si cell
SAMPLE_QE_OUTPUT = """
     Program PWSCF v.6.5 starts on 10Jan2024 at 10:00:00

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
     "P. Giannozzi et al., J. Phys.: Condens. Matter 21 395502 (2009);
     URL http://www.quantum-espresso.org",
     in publications or presentations arising from this work.

     Parallel version (MPI), running on     4 processors

     ... [various setup info] ...

     Self-consistent field calculation ...

     ... [convergence info] ...

!    total energy              =     -15.85217439 Ry

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =     -0.00000014    -0.00000014    -0.00000014
     atom    2   force =      0.00000014     0.00000014     0.00000014

     Total force =     0.00000000     Total SCF correction =     0.00000000

     ... [other info] ...

     total   stress  (Ry/bohr**3)                (kbar)     P=      -0.34
      -0.00000215    -0.00000000     0.00000000
      -0.00000000    -0.00000215     0.00000000
       0.00000000     0.00000000    -0.00000215

     JOB DONE.
"""


def test_parse_qe_output_success():
    """Tests successful parsing of energy, forces, and stress from QE output."""
    results = parse_qe_output(SAMPLE_QE_OUTPUT)

    assert "energy" in results
    assert "forces" in results
    assert "stress" in results

    # Check energy (converted from Ry to eV)
    expected_energy = -15.85217439 * 13.605693122994
    assert results["energy"] == pytest.approx(expected_energy, rel=1e-6)

    # Check forces (converted from Ry/au to eV/A)
    expected_forces_ry_au = np.array(
        [
            [-0.00000014, -0.00000014, -0.00000014],
            [0.00000014, 0.00000014, 0.00000014],
        ]
    )
    ry_au_to_ev_a = 13.605693122994 / 0.529177210903
    expected_forces_ev_a = expected_forces_ry_au * ry_au_to_ev_a
    assert np.allclose(results["forces"], expected_forces_ev_a, rtol=1e-6)
    assert results["forces"].shape == (2, 3)

    # Check stress (converted from kbar to eV/A^3)
    expected_stress_kbar = np.array(
        [
            [-0.00000215, -0.00000000, 0.00000000],
            [-0.00000000, -0.00000215, 0.00000000],
            [0.00000000, 0.00000000, -0.00000215],
        ]
    )
    kbar_to_ev_a3 = 1 / 160.21766208
    expected_stress_ev_a3 = expected_stress_kbar * kbar_to_ev_a3
    assert np.allclose(results["stress"], expected_stress_ev_a3, rtol=1e-6)
    assert results["stress"].shape == (3, 3)


def test_parse_qe_output_failure():
    """Tests that parsing returns an empty dict for invalid output."""
    invalid_output = "This is not a Quantum Espresso output file."
    results = parse_qe_output(invalid_output)
    assert results == {}
