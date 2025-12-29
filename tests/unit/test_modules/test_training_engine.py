from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db():
    """
    Fixture for a mocked AseDB that returns pre-defined data, including an
    Atoms object with a calculator attached.
    """
    db = MagicMock(spec=AseDB)
    atoms = Atoms('H', positions=[[0, 0, 0]])
    dft_energy = -10.0
    dft_forces = np.array([[0.1, 0.2, 0.3]])

    # Attach a calculator, just as the LabellingEngine would
    calc = SinglePointCalculator(atoms, energy=dft_energy, forces=dft_forces)
    atoms.calc = calc

    dft_kvp = { "was_successful": True }
    db.get.return_value = (atoms, dft_kvp)
    return db

@pytest.fixture
def training_config():
    """Fixture for a standard TrainingConfig."""
    return TrainingConfig(
        model_type='mace',
        learning_rate=0.01,
        num_epochs=10,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential='lj'
    )

@patch('mlip_autopipec.utils.baseline_potentials.calculate_lj')
def test_training_engine_load_and_prepare_data_delta_learning(
    mock_calculate_lj, mock_db, training_config
):
    """
    Tests the data loading and preparation step with delta learning enabled.
    Verifies that the baseline potential is called and that the resulting
    training data contains the correct DFT-baseline residual.
    """
    # Arrange
    baseline_energy = -8.0
    baseline_forces = np.array([[0.05, 0.1, 0.15]])
    mock_calculate_lj.return_value = (baseline_energy, baseline_forces)

    engine = TrainingEngine(config=training_config, db=mock_db)

    # Act
    prepared_data = engine._load_and_prepare_data(ids=[1])

    # Assert
    mock_db.get.assert_called_once_with(1)
    mock_calculate_lj.assert_called_once()

    atomic_data = prepared_data[0]
    expected_energy_delta = -10.0 - baseline_energy
    expected_forces_delta = np.array([[0.1, 0.2, 0.3]]) - baseline_forces

    assert np.isclose(atomic_data.energy.item(), expected_energy_delta)
    assert np.allclose(atomic_data.forces.numpy(), expected_forces_delta)
