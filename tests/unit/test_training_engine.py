import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from ase.build import bulk
from src.mlip_autopipec.modules.d_training_engine import TrainingEngine
from src.mlip_autopipec.data.models import TrainingConfig
from src.mlip_autopipec.data.database import AseDB

@pytest.fixture
def mock_db_with_data(tmp_path):
    """Fixture for a mocked AseDB that returns a sample data point."""
    db_path = tmp_path / "test.db"
    db = AseDB(str(db_path))

    atoms = bulk("Si", "diamond", a=5.43)
    dft_forces = np.random.rand(len(atoms), 3)
    dft_energy = -100.0

    atoms.calc = MagicMock()
    atoms.get_potential_energy = MagicMock(return_value=dft_energy)
    atoms.get_forces = MagicMock(return_value=dft_forces)

    db.get = MagicMock(return_value=atoms)
    return db, dft_energy, dft_forces

@patch("src.mlip_autopipec.utils.baseline_potentials.lennard_jones_potential")
def test_training_engine_data_preparation(mock_lj, mock_db_with_data):
    """
    Tests that the TrainingEngine's internal data preparation method correctly
    calculates the delta energy and forces for Delta Learning.
    """
    db, dft_energy, dft_forces = mock_db_with_data

    baseline_energy = -80.0
    baseline_forces = np.ones((2, 3)) * 0.1
    mock_lj.return_value = (baseline_energy, baseline_forces)

    config = TrainingConfig(
        model_type="mace",
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lennard_jones",
    )

    engine = TrainingEngine(config=config, db=db, epsilon_lj=0.1, sigma_lj=3.0)

    # Directly test the private method
    prepared_atoms_list = engine._load_and_prepare_data(ids=[1])

    db.get.assert_called_with(1)
    mock_lj.assert_called_once()

    # Verify the delta energy and forces stored in the atoms.info dict
    processed_atoms = prepared_atoms_list[0]
    delta_energy = processed_atoms.info['energy']
    delta_forces = processed_atoms.info['forces']

    expected_delta_energy = dft_energy - baseline_energy
    expected_delta_forces = dft_forces - baseline_forces

    assert delta_energy == pytest.approx(expected_delta_energy)
    assert np.allclose(delta_forces, expected_delta_forces)
