import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from ase.build import bulk
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig, DFTResult
from mlip_autopipec.modules.d_training_engine import TrainingEngine

@pytest.fixture
def mock_db_with_data(tmp_path):
    """A fixture to create a mock AseDB with some test data."""
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)

    atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        total_energy_ev=-150.0,
        forces=np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]]),
        stress=np.zeros((3,3)),
        was_successful=True,
    )
    db.write(atoms, dft_result)
    return db

def test_load_and_prepare_data_delta_learning(mock_db_with_data):
    """
    Tests that the data loading correctly prepares data for delta learning.
    """
    config = TrainingConfig(
        model_type="mace",
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lennard_jones",
    )
    engine = TrainingEngine(config, mock_db_with_data)

    prepared_data = engine._load_and_prepare_data([1])

    # 1. Check that we have one data point
    assert len(prepared_data) == 1

    # 2. Check the structure of the data point
    data_point = prepared_data[0]
    assert "atoms" in data_point
    assert "energy_delta" in data_point
    assert "forces_delta" in data_point

    # 3. Verify the delta calculation
    # In this case, the baseline LJ potential for Si will be non-zero
    atoms = mock_db_with_data.get_atoms(1)
    from mlip_autopipec.utils.baseline_potentials import lennard_jones_potential
    base_energy, base_forces = lennard_jones_potential(atoms)

    expected_energy_delta = -150.0 - base_energy
    expected_forces_delta = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]]) - base_forces

    assert data_point["energy_delta"] == pytest.approx(expected_energy_delta)
    assert np.allclose(data_point["forces_delta"], expected_forces_delta)

@patch("mlip_autopipec.modules.d_training_engine.run_train_cli")
def test_training_engine_execute_call(mock_run_train_cli, mock_db_with_data):
    """
    Tests that the execute method calls the underlying MACE training function
    with the correct parameters.
    """
    config = TrainingConfig(
        model_type="mace",
        learning_rate=0.01,
        num_epochs=10,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lennard_jones",
    )
    engine = TrainingEngine(config, mock_db_with_data)

    model_path = engine.execute([1])

    assert "mace_training_final.pt" in model_path
    mock_run_train_cli.assert_called_once()

    # Check that some key arguments were passed to the MACE trainer
    args, kwargs = mock_run_train_cli.call_args
    train_args = args[0]
    assert train_args.name == "mace_training"
    assert train_args.train_file is not None # Path to a temporary file
    assert train_args.learning_rate == 0.01
    assert train_args.num_epochs == 10
    assert "'energy_key': 'energy_delta'" in train_args.config_type_weights
    assert "'forces_key': 'forces_delta'" in train_args.config_type_weights
