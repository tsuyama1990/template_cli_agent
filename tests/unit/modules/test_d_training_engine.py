"""Unit tests for the TrainingEngine module."""
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from pathlib import Path
from ase import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig, DFTResult
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db(tmp_path: Path) -> MagicMock:
    """Provides a mock AseDB instance."""
    db = MagicMock(spec=AseDB)

    # Setup mock data for a successful calculation
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    dft_result = DFTResult(
        total_energy_ev=-31.0,
        forces=[[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]],
        stress=[[0.0]*3]*3,
        was_successful=True
    )
    db.get.return_value = (atoms, dft_result)

    return db


@pytest.fixture
def training_config() -> TrainingConfig:
    """Provides a standard training configuration."""
    return TrainingConfig(
        model_type="mace",
        learning_rate=1e-3,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lj"
    )

@patch("mlip_autopipec.modules.d_training_engine.baseline_potentials.calculate_lennard_jones")
def test_load_and_prepare_data_delta_learning(mock_calculate_lj, mock_db: MagicMock, training_config: TrainingConfig):
    """
    Tests that the data preparation step correctly calculates the delta
    between the DFT results and the baseline potential.
    """
    # Arrange
    # Mock the baseline potential to return predictable values
    baseline_energy = -1.0
    baseline_forces = np.array([[0.0, 0.0, 0.05], [0.0, 0.0, -0.05]])
    mock_calculate_lj.return_value = (baseline_energy, baseline_forces)

    engine = TrainingEngine(config=training_config, db=mock_db)

    # Act
    prepared_data = engine._load_and_prepare_data(ids=[1])

    # Assert
    assert len(prepared_data) == 1

    # Check that the baseline potential was called
    mock_calculate_lj.assert_called_once()

    # Check that the target energy and forces are the *delta* values
    atoms_out = prepared_data[0]
    expected_delta_energy = -31.0 - baseline_energy  # -30.0
    expected_delta_forces = np.array([[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]]) - baseline_forces

    assert atoms_out.info['energy'] == pytest.approx(expected_delta_energy)
    assert np.allclose(atoms_out.info['forces'], expected_delta_forces)


@patch("mlip_autopipec.modules.d_training_engine.baseline_potentials.calculate_lennard_jones")
def test_load_and_prepare_data_no_delta_learning(mock_calculate_lj, mock_db: MagicMock, training_config: TrainingConfig):
    """
    Tests that the data preparation step uses raw DFT values when delta learning is disabled.
    """
    # Arrange
    training_config.delta_learn = False
    engine = TrainingEngine(config=training_config, db=mock_db)

    # Act
    prepared_data = engine._load_and_prepare_data(ids=[1])

    # Assert
    assert len(prepared_data) == 1

    # Check that the baseline potential was NOT called
    mock_calculate_lj.assert_not_called()

    # Check that the target energy and forces are the raw DFT values
    atoms_out = prepared_data[0]
    raw_dft_result = mock_db.get.return_value[1]

    assert atoms_out.info['energy'] == pytest.approx(raw_dft_result.total_energy_ev)
    assert np.allclose(atoms_out.info['forces'], raw_dft_result.forces)
