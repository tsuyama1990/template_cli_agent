"""Unit tests for the LabellingEngine."""

from unittest.mock import MagicMock, Mock

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.modules.labelling_engine import LabellingEngine


@pytest.fixture
def mock_db():
    """Fixture for a mocked AseDB."""
    return MagicMock()


@pytest.fixture
def mock_runner():
    """Fixture for a mocked ProcessRunner."""
    return MagicMock()


@pytest.fixture
def config():
    """Fixture for a sample configuration dictionary."""
    return {
        "dft_command": "pw.x -in stdin -out stdout",
        "pseudo_dir": "/path/to/pseudos",
        "ecutwfc": 60.0,
    }


# A more realistic sample output that matches the parser's regex
SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =     -11.45567302 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000    0.00000000
     atom    2   force =     0.00000101    0.00000000    0.00000000

total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
     0.10000000   0.40000000   0.50000000
     0.40000000   0.20000000   0.60000000
     0.50000000   0.60000000   0.30000000
"""

SAMPLE_FAILURE_OUTPUT = "convergence NOT achieved"


def test_labelling_engine_run_success(mock_db, mock_runner, config):
    """Test a successful run of the LabellingEngine."""
    # Arrange
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], info={"db_id": 1})
    mock_db.get_atoms_by_state.return_value = [atoms]
    mock_runner.run.return_value = SAMPLE_SUCCESS_OUTPUT

    engine = LabellingEngine(mock_db, mock_runner, config)

    # Act
    engine.run()

    # Assert
    mock_db.get_atoms_by_state.assert_called_once_with("unlabelled")
    mock_runner.run.assert_called_once()
    # Check that update_with_dft_results was called with a DFTResult object
    mock_db.update_with_dft_results.assert_called_once()
    call_args = mock_db.update_with_dft_results.call_args[0]
    assert call_args[0] == 1  # atoms_id
    assert isinstance(call_args[1], DFTResult)
    assert call_args[1].energy == -11.45567302


def test_labelling_engine_run_dft_failure(mock_db, mock_runner, config):
    """Test a run where the DFT calculation fails."""
    # Arrange
    atoms = Atoms("H", info={"db_id": 1})
    mock_db.get_atoms_by_state.return_value = [atoms]
    mock_runner.run.return_value = SAMPLE_FAILURE_OUTPUT

    engine = LabellingEngine(mock_db, mock_runner, config)

    # Act
    engine.run()

    # Assert
    mock_db.get_atoms_by_state.assert_called_once_with("unlabelled")
    mock_runner.run.assert_called_once()
    # Crucially, assert that the database was NOT updated with bad data
    mock_db.update_with_dft_results.assert_not_called()


def test_labelling_engine_handles_empty_database(mock_db, mock_runner, config):
    """Test that the engine does nothing if there are no unlabelled atoms."""
    # Arrange
    mock_db.get_atoms_by_state.return_value = []
    engine = LabellingEngine(mock_db, mock_runner, config)

    # Act
    engine.run()

    # Assert
    mock_db.get_atoms_by_state.assert_called_once_with("unlabelled")
    mock_runner.run.assert_not_called()
    mock_db.update_with_dft_results.assert_not_called()


def test_generate_qe_input_content(config):
    """Test the content of the generated QE input."""
    engine = LabellingEngine(Mock(), Mock(), config)
    atoms = Atoms("Si", cell=np.eye(3))
    input_str = engine._generate_qe_input(atoms)

    # Make assertions less strict to whitespace and quotes
    assert "calculation" in input_str
    assert "'scf'" in input_str
    assert "ecutwfc" in input_str
    assert str(config["ecutwfc"]) in input_str
    assert "nat" in input_str
    assert "= 1" in input_str
    assert "ATOMIC_SPECIES" in input_str
    assert "Si" in input_str
    assert "ATOMIC_POSITIONS" in input_str
    assert "CELL_PARAMETERS" in input_str
