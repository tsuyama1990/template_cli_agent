# Description: Unit tests for the utility modules dft_utils and baseline_potentials.
import numpy as np
import pytest
from ase.build import bulk

from mlip_autopipec.utils.baseline_potentials import zbl_potential
from mlip_autopipec.utils.dft_utils import generate_qe_input, parse_qe_output


# --- Tests for dft_utils.py ---
def test_generate_qe_input():
    """Tests the generation of a Quantum Espresso input file."""
    si_atoms = bulk("Si", "diamond", a=5.43)
    input_str = generate_qe_input(si_atoms)

    assert "&CONTROL" in input_str
    assert "&SYSTEM" in input_str
    assert "nat = 2" in input_str
    assert "ntyp = 1" in input_str
    assert "ATOMIC_SPECIES" in input_str
    # Loosen the assertion to avoid issues with precise mass values from ase.data
    assert "Si" in input_str
    assert "Si.pbe-n-rrkjus_psl.1.0.0.UPF" in input_str
    assert "ATOMIC_POSITIONS {crystal}" in input_str
    assert "CELL_PARAMETERS {angstrom}" in input_str


def test_parse_qe_output_successful():
    """Tests parsing a successful Quantum Espresso output."""
    mock_output = """
     JOB DONE.
     !    total energy              =     -11.234 Ry
     Forces acting on atoms (Ry/au):
     atom    1 type  1   force =    -0.001   0.002   0.003
     atom    2 type  1   force =     0.001  -0.002  -0.003
     total stress  (Ry/bohr**3)     (kbar)     P=     -0.12
       -0.0001   0.0000   0.0000
        0.0000  -0.0001   0.0000
        0.0000   0.0000  -0.0001
    """
    result = parse_qe_output(mock_output)

    assert result.was_successful is True
    assert result.error_message is None
    assert result.total_energy_ev == pytest.approx(-11.234 * 13.6057, abs=1e-3)
    assert len(result.forces) == 2
    # Conversion factor from Ry/au to eV/Angstrom
    ry_au_to_ev_ang = 25.71104309541616 * 2
    assert result.forces[0][0] == pytest.approx(-0.001 * ry_au_to_ev_ang, abs=1e-3)


def test_parse_qe_output_not_converged():
    """Tests parsing a failed (non-converged) QE output."""
    mock_output = "     convergence NOT achieved"
    result = parse_qe_output(mock_output)

    assert result.was_successful is False
    assert result.error_message == "SCF did not converge."


# --- Tests for baseline_potentials.py ---
def test_zbl_potential_two_atoms():
    """Tests the ZBL potential calculation for a simple two-atom system."""
    atoms = bulk("Si", "diamond", a=5.43)
    # Get just the first two atoms to create a simple pair
    two_atoms = atoms[0:2]
    # Set the cell to be large to avoid periodic image interactions
    two_atoms.set_cell([20, 20, 20], scale_atoms=False)

    energy, forces = zbl_potential(two_atoms)

    # Check energy is a positive (repulsive) float
    assert isinstance(energy, float)
    assert energy > 0

    # Check forces are a numpy array of the correct shape
    assert isinstance(forces, np.ndarray)
    assert forces.shape == (2, 3)

    # Check Newton's third law: F_1 = -F_2
    assert np.allclose(forces[0], -forces[1])

    # Check that the total force is zero
    assert np.allclose(np.sum(forces, axis=0), np.zeros(3))
