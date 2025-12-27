import pytest
import numpy as np
import torch
from unittest.mock import patch, MagicMock, ANY
from ase import Atoms
import unittest

from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult, TrainingConfig

@pytest.fixture
def mock_db_with_data(tmp_path):
    """Fixture for a mock AseDB with one successful entry."""
    db_path = tmp_path / "test.db"
    db = AseDB(str(db_path))
    atoms = Atoms('H', positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    # Use a real calculator to store data correctly
    from ase.calculators.singlepoint import SinglePointCalculator
    calc = SinglePointCalculator(atoms, energy=-10.0, forces=np.array([[0.1, 0.2, 0.3]]))
    atoms.calc = calc
    db_id = db.write(atoms, DFTResult(was_successful=True, total_energy_ev=-10.0, forces=[[0.1,0.2,0.3]], stress=[[0.,0.,0.]]*3))
    return db, db_id

@pytest.fixture
def training_config():
    """Fixture for a standard TrainingConfig."""
    return TrainingConfig(
        model_type='mace',
        learning_rate=1e-3,
        num_epochs=2,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential='lj'
    )

def test_load_and_prepare_data_delta_learn(mock_db_with_data, training_config):
    """Tests that data preparation correctly applies the delta."""
    db, db_id = mock_db_with_data
    engine = TrainingEngine(config=training_config, db=db) # delta_learn=True

    prepared_atoms_list = engine._load_and_prepare_data([db_id])

    assert len(prepared_atoms_list) == 1
    atoms = prepared_atoms_list[0]

    # LJ for a single atom is zero energy and forces
    expected_energy = -10.0
    expected_forces = np.array([[0.1, 0.2, 0.3]])

    assert np.isclose(atoms.get_potential_energy(), expected_energy)
    assert np.allclose(atoms.get_forces(), expected_forces)

def test_load_and_prepare_data_no_delta(mock_db_with_data, training_config):
    """Tests that data preparation passes DFT values through when delta is false."""
    db, db_id = mock_db_with_data
    training_config.delta_learn = False
    engine = TrainingEngine(config=training_config, db=db)

    prepared_atoms_list = engine._load_and_prepare_data([db_id])
    atoms = prepared_atoms_list[0]

    assert np.isclose(atoms.get_potential_energy(), -10.0)
    assert np.allclose(atoms.get_forces(), np.array([[0.1, 0.2, 0.3]]))


@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch.optim.Adam')
@patch('mlip_autopipec.modules.d_training_engine.torch.save')
@patch('mlip_autopipec.modules.d_training_engine.WeightedEnergyForcesLoss')
def test_training_engine_execute(mock_loss_fn, mock_torch_save, mock_adam, mock_mace, mock_db_with_data, training_config):
    """Tests the main execute method, mocking all ML components."""
    db, db_id = mock_db_with_data

    # Configure mocks
    mock_model = MagicMock()
    mock_mace.return_value = mock_model
    mock_model.return_value = {'energy': torch.tensor([1.0]), 'forces': torch.tensor([[1.0, 1.0, 1.0]])}

    mock_loss = MagicMock()
    mock_loss_fn.return_value = mock_loss

    # Create a mock tensor with a mock backward method and item method
    mock_loss_tensor = MagicMock(spec=torch.Tensor)
    mock_loss_tensor.backward = MagicMock()
    mock_loss_tensor.item.return_value = 0.5
    mock_loss.return_value = mock_loss_tensor

    engine = TrainingEngine(config=training_config, db=db)
    model_path = engine.execute([db_id])

    # Assertions
    mock_mace.assert_called_once()
    mock_adam.assert_called_once_with(ANY, lr=training_config.learning_rate)

    # Check that the training loop ran
    assert mock_loss.call_count == training_config.num_epochs
    assert mock_loss_tensor.backward.call_count == training_config.num_epochs

    mock_torch_save.assert_called_once_with(ANY, model_path)
    assert model_path.startswith("models/")
    assert model_path.endswith(".pt")
