"""Unit tests for the Training Engine."""
import pytest
from unittest.mock import patch, MagicMock
import numpy as np
from ase.atoms import Atoms
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.data.database import AseDB

@pytest.fixture
def mock_db():
    db = MagicMock(spec=AseDB)
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
    db_entry = {
        "atoms": atoms,
        "total_energy_ev": -10.0,
        "forces": np.array([[0, 0, 0.1], [0, 0, -0.1]]),
        "was_successful": True,
    }
    db.get.return_value = db_entry
    return db

@pytest.fixture
def training_config():
    return TrainingConfig(
        model_type='mace',
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential='lennard_jones'
    )

def mock_tensor_creation(data, dtype):
    """Custom mock for torch.tensor."""
    if isinstance(data, (float, int)):
        return data  # Return raw value for energy
    else:
        # Return a mock tensor with a .numpy() method for forces
        mock_tensor = MagicMock()
        mock_tensor.numpy.return_value = np.array(data)
        return mock_tensor

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch')
@patch('mlip_autopipec.modules.d_training_engine.Adam')
@patch('mlip_autopipec.modules.d_training_engine.AtomicNumberTable')
@patch('mlip_autopipec.modules.d_training_engine.config_from_atoms')
@patch('mlip_autopipec.modules.d_training_engine.AtomicData')
def test_training_engine_delta_learning(mock_atomic_data, mock_config, mock_z_table, mock_adam, mock_torch, mock_mace, mock_db, training_config):
    """Test the data preparation with delta learning enabled."""
    mock_torch.tensor.side_effect = mock_tensor_creation

    class MockAtomicData:
        def to_device(self, device): return self
    mock_atomic_data.from_config.return_value = MockAtomicData()

    engine = TrainingEngine(config=training_config, db=mock_db)
    prepared_data = engine._load_and_prepare_data(ids=[1])

    db_entry = mock_db.get()
    dft_energy = db_entry["total_energy_ev"]
    dft_forces = db_entry["forces"]

    assert prepared_data[0].energy < dft_energy
    assert not np.allclose(prepared_data[0].forces.numpy(), dft_forces)

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch')
@patch('mlip_autopipec.modules.d_training_engine.Adam')
@patch('mlip_autopipec.modules.d_training_engine.AtomicNumberTable')
@patch('mlip_autopipec.modules.d_training_engine.config_from_atoms')
@patch('mlip_autopipec.modules.d_training_engine.AtomicData')
def test_training_engine_no_delta_learning(mock_atomic_data, mock_config, mock_z_table, mock_adam, mock_torch, mock_mace, mock_db, training_config):
    """Test the data preparation with delta learning disabled."""
    training_config.delta_learn = False
    mock_torch.tensor.side_effect = mock_tensor_creation

    class MockAtomicData:
        def to_device(self, device): return self
    mock_atomic_data.from_config.return_value = MockAtomicData()

    engine = TrainingEngine(config=training_config, db=mock_db)
    prepared_data = engine._load_and_prepare_data(ids=[1])

    db_entry = mock_db.get()
    dft_energy = db_entry["total_energy_ev"]
    dft_forces = db_entry["forces"]

    assert np.isclose(prepared_data[0].energy, dft_energy)
    assert np.allclose(prepared_data[0].forces.numpy(), dft_forces)
