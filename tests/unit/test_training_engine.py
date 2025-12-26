import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig

@pytest.fixture
def mock_db_with_data(tmp_path):
    db_path = tmp_path / "test.db"
    db = AseDB(str(db_path))
    atoms = bulk("Si", "diamond", a=5.43)

    # Attach results via a calculator
    calc = SinglePointCalculator(
        atoms,
        energy=-100.0,
        forces=np.array([[0.1, 0.1, 0.1]] * 2),
        stress=np.zeros((3, 3))
    )
    atoms.calc = calc

    db.write(atoms, metadata={'was_successful': True})
    return db, -100.0

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch')
def test_training_engine_delta_learn(mock_torch, mock_mace, mock_db_with_data):
    mock_db, dft_energy = mock_db_with_data
    training_config = TrainingConfig(
        model_type='mace',
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential='zbl'
    )

    engine = TrainingEngine(config=training_config, db=mock_db)
    prepared_data = engine._load_and_prepare_data(ids=[1])

    assert prepared_data[0].energy != dft_energy

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch')
def test_training_engine_direct_learn(mock_torch, mock_mace, mock_db_with_data):
    mock_db, dft_energy = mock_db_with_data
    training_config = TrainingConfig(
        model_type='mace',
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=False,
        baseline_potential=''
    )

    engine = TrainingEngine(config=training_config, db=mock_db)
    prepared_data = engine._load_and_prepare_data(ids=[1])

    assert prepared_data[0].energy == dft_energy
