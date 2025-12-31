import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.config.models import DFTParams
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine, QEOutputParsingError


@pytest.fixture
def sample_dft_params():
    """Provides a sample DFTParams object for tests."""
    return DFTParams(
        command="pw.x",
        pseudopotentials={"Si": "Si.upf"},
        ecutwfc=60.0,
        kpoints_density=5000,
    )


@pytest.fixture
def sample_atoms():
    """Provides a sample Atoms object for a silicon dimer."""
    return Atoms("Si2", positions=[[0, 0, 0], [2.35, 0, 0]], cell=[10, 10, 10], pbc=True)


def test_generate_qe_input(sample_dft_params, sample_atoms):
    """
    Tests that a valid Quantum Espresso input string is generated robustly.
    """
    engine = LabellingEngine(dft_config=sample_dft_params)
    input_str = engine._generate_qe_input(sample_atoms)
    lines = [line.strip() for line in input_str.split('\n')]

    def find_line_starting_with(prefix):
        """Helper to find the index of the first line starting with a prefix."""
        for i, line in enumerate(lines):
            if line.startswith(prefix):
                return i
        return -1

    # Check for presence of key cards using the robust helper
    assert "&CONTROL" in lines
    assert "&SYSTEM" in lines
    assert "&ELECTRONS" in lines
    assert "ATOMIC_SPECIES" in lines
    assert find_line_starting_with("ATOMIC_POSITIONS") != -1
    assert find_line_starting_with("CELL_PARAMETERS") != -1
    assert find_line_starting_with("K_POINTS") != -1

    # Verify a specific parameter to ensure config is used
    found_ecut = any("ecutwfc" in line and "60.0" in line for line in lines)
    assert found_ecut, "ecutwfc not found or incorrect in output."

    # Verify the atomic species line robustly
    species_index = lines.index("ATOMIC_SPECIES")
    species_line = lines[species_index + 1]
    assert species_line.startswith("Si")
    assert species_line.endswith("Si.upf")

    # Verify one of the atomic positions robustly
    pos_index = find_line_starting_with("ATOMIC_POSITIONS")
    pos_line = lines[pos_index + 2]  # Second atom
    assert pos_line.startswith("Si")
    coords = [float(x) for x in pos_line.split()[1:]]
    assert coords == pytest.approx([2.35, 0.0, 0.0])


def load_qe_output(filename: str) -> str:
    """Helper to load a QE output file from the test data directory."""
    path = Path(__file__).parent / ".." / "test_data" / "qe_outputs" / filename
    return path.read_text()


def test_parse_qe_output_success(sample_dft_params):
    """
    Tests successful parsing of a complete and valid QE output file.
    """
    engine = LabellingEngine(dft_config=sample_dft_params)
    output_text = load_qe_output("successful_run.out")
    results = engine._parse_qe_output(output_text)

    # Check energy (QE output is in Ry, should be converted to eV)
    assert results["energy"] == pytest.approx(-15.85311181 * 13.605693122994)

    # Check forces (QE output is in Ry/au, should be converted to eV/Ang)
    expected_forces = np.array([
        [-0.00123456, 0.0, 0.0],
        [0.00123456, 0.0, 0.0]
    ]) * (13.605693122994 / 0.529177210903)
    assert np.allclose(results["forces"], expected_forces)

    # Check stress (QE output is in kbar, should be converted to eV/Ang^3)
    # Voigt order: xx, yy, zz, yz, xz, xy
    expected_stress_kbar = np.array([-0.20845347, -0.17647228, -0.17647228, 0.0, 0.0, 0.0])
    KBAR_TO_EV_ANG3 = 1.0 / 160.21766208 * 0.1
    expected_stress_ev_ang3 = expected_stress_kbar * KBAR_TO_EV_ANG3
    assert np.allclose(results["stress"], expected_stress_ev_ang3)


def test_parse_qe_output_no_energy(sample_dft_params):
    """Tests that an error is raised if the final energy is missing."""
    engine = LabellingEngine(dft_config=sample_dft_params)
    output_text = "Some output without the final energy."
    with pytest.raises(QEOutputParsingError, match="Final energy not found"):
        engine._parse_qe_output(output_text)


def test_parse_qe_output_no_forces(sample_dft_params):
    """Tests that an error is raised if the forces block is missing."""
    engine = LabellingEngine(dft_config=sample_dft_params)
    output_text = "!    total energy              =     -15.85311181 Ry"
    with pytest.raises(QEOutputParsingError, match="Forces block not found"):
        engine._parse_qe_output(output_text)


def test_parse_qe_output_no_stress(sample_dft_params):
    """Tests that an error is raised if the stress block is missing."""
    engine = LabellingEngine(dft_config=sample_dft_params)
    output_text = """
    !    total energy              =     -15.85311181 Ry
    Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =    -0.00123456   0.00000000   0.00000000
     atom    2 type  1   force =     0.00123456   0.00000000   0.00000000

    """
    with pytest.raises(QEOutputParsingError, match="Stress block not found"):
        engine._parse_qe_output(output_text)


@patch("subprocess.run")
def test_execute_qe_success(mock_subprocess_run, sample_dft_params):
    """Tests the _execute_qe method for a successful execution."""
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.stdout = "Successful QE run output."
    mock_process.stderr = ""
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(dft_config=sample_dft_params)
    success, stdout, stderr = engine._execute_qe("fake_input")

    assert success is True
    assert stdout == "Successful QE run output."
    assert stderr == ""
    mock_subprocess_run.assert_called_once_with(
        sample_dft_params.command.split(),
        input="fake_input",
        text=True,
        capture_output=True,
        check=False,
        timeout=300,
    )


@patch("subprocess.run")
def test_execute_qe_failure(mock_subprocess_run, sample_dft_params):
    """Tests the _execute_qe method for a failed execution."""
    mock_process = MagicMock()
    mock_process.returncode = 1
    mock_process.stdout = ""
    mock_process.stderr = "QE calculation failed."
    mock_subprocess_run.return_value = mock_process

    engine = LabellingEngine(dft_config=sample_dft_params)
    success, stdout, stderr = engine._execute_qe("fake_input")

    assert success is False
    assert stderr == "QE calculation failed."


@patch("subprocess.run", side_effect=FileNotFoundError)
def test_execute_qe_command_not_found(mock_subprocess_run, sample_dft_params):
    """Tests handling of FileNotFoundError."""
    engine = LabellingEngine(dft_config=sample_dft_params)
    success, stdout, stderr = engine._execute_qe("fake_input")

    assert success is False
    assert "Command 'pw.x' not found" in stderr


@patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="pw.x", timeout=300))
def test_execute_qe_timeout(mock_subprocess_run, sample_dft_params):
    """Tests handling of TimeoutExpired."""
    engine = LabellingEngine(dft_config=sample_dft_params)
    success, stdout, stderr = engine._execute_qe("fake_input")

    assert success is False
    assert "Quantum Espresso command timed out" in stderr

# More detailed tests will be added as the implementation proceeds.
# For now, we are just setting up the test structure.
