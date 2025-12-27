import pytest
import numpy as np
from ase.build import bulk
from src.mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output

# Sample successful Quantum Espresso output
SAMPLE_SUCCESS_OUTPUT = """
     crystal axes: (cart. coord. in Angstrom)
          a(1) = (   0.000000   2.715000   2.715000  )
          a(2) = (   2.715000   0.000000   2.715000  )
          a(3) = (   2.715000   2.715000   0.000000  )

     site n.     atom                  positions (cryst. coord.)
        1           Si  tau(   1) = (   0.000000000   0.000000000   0.000000000  )
        2           Si  tau(   2) = (   0.250000000   0.250000000   0.250000000  )

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00000002    -0.00000002    -0.00000002
     atom    2 type  1   force =     0.00000002     0.00000002     0.00000002

     Total Force =     0.00000004     Total SCF correction =     0.00000000


!    total energy              =    -111.92138384 Ry

          total stress  (Ry/bohr**3)            (kbar)       P=   -2.29
   -0.00001452   0.00000000   0.00000000      -0.23      0.00      0.00
    0.00000000  -0.00001452   0.00000000       0.00     -0.23      0.00
    0.00000000   0.00000000  -0.00001452       0.00      0.00     -0.23

"""

# Sample failed (non-converged) Quantum Espresso output
SAMPLE_FAILURE_OUTPUT = """
     Maximum number of SCF cycles reached.
     Message from routine cdiaghg:
     not converged in 100 iterations

     End of self-consistent calculation
"""

def test_generate_qe_input():
    """Tests the generation of a Quantum Espresso input file string."""
    atoms = bulk("Si", "diamond", a=5.43)
    pseudo_dir = "/path/to/pseudos"
    pseudopotentials = {"Si": "Si.upf"}
    kpoints = (4, 4, 4)
    ecutwfc = 60.0

    input_str = generate_qe_input(
        atoms=atoms,
        pseudo_dir=pseudo_dir,
        pseudopotentials=pseudopotentials,
        kpoints=kpoints,
        ecutwfc=ecutwfc,
    )

    assert "ATOMIC_SPECIES" in input_str
    assert "Si  28.0850  Si.upf" in input_str
    assert "ATOMIC_POSITIONS" in input_str
    assert "CELL_PARAMETERS" in input_str
    assert "K_POINTS" in input_str
    assert "4 4 4 0 0 0" in input_str
    assert "ecutwfc = 60.0" in input_str
    assert "nat = 2" in input_str

def test_parse_qe_output_success():
    """Tests parsing a successful QE output."""
    result = parse_qe_output(SAMPLE_SUCCESS_OUTPUT)

    assert result.was_successful
    assert result.error_message is None
    assert result.total_energy_ev == pytest.approx(-111.92138384 * 13.6057) # Ry to eV

    # Check forces
    assert len(result.forces) == 2
    # Conversion from Ry/au to eV/Angstrom
    ry_au_to_ev_a = 13.6057 / 0.529177
    assert np.allclose(result.forces[0], [-0.00000002 * ry_au_to_ev_a] * 3)

    # Check stress
    # Conversion from kbar to eV/Angstrom^3
    kbar_to_ev_a3 = 1.0 / 160.21766208
    expected_stress_diag = -0.23 * kbar_to_ev_a3
    expected_stress = np.diag([expected_stress_diag] * 3)
    assert np.allclose(result.stress, expected_stress)


def test_parse_qe_output_failure():
    """Tests parsing a failed (non-converged) QE output."""
    result = parse_qe_output(SAMPLE_FAILURE_OUTPUT)

    assert not result.was_successful
    assert "SCF did not converge" in result.error_message
    assert result.total_energy_ev == 0.0
    assert result.forces == []
    assert result.stress == []
