import pytest
from unittest.mock import MagicMock
from ase.build import bulk
import numpy as np

from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential

# Mock AtomsRow object to simulate what comes out of ase.db
class MockAtomsRow:
    def __init__(self, atoms, data, key_value_pairs):
        self._atoms = atoms
        self.data = data
        self.key_value_pairs = key_value_pairs

    def toatoms(self):
        return self._atoms

@pytest.fixture
def mock_db():
    """Fixture for a mocked AseDB instance with a sample data point."""
    db = MagicMock()

    atoms = bulk("Ar", "fcc", a=5.2)
    dft_data = {
        "total_energy_ev": -1.0,
        "forces": [[0.1, 0.0, 0.0]] * 4,
        "stress": [],
        "was_successful": True,
        "error_message": None
    }
    kvp = {"was_successful": True, "total_energy_ev": -1.0}

    mock_row = MockAtomsRow(atoms, dft_data, kvp)
    db.get.return_value = mock_row
    return db

def test_training_engine_load_data_no_delta(mock_db):
    """
    Tests the data loading logic when Delta Learning is disabled.
    """
    config = TrainingConfig(
        model_type="mace", learning_rate=0.01, num_epochs=1,
        r_cut=5.0, delta_learn=False, baseline_potential="lj"
    )
    engine = TrainingEngine(config, mock_db)

    prepared_data = engine._load_and_prepare_data(ids=[1])

    assert len(prepared_data) == 1

    data_point = prepared_data[0]
    # Without delta learning, targets should be the raw DFT values
    assert data_point["target_energy"] == -1.0
    assert np.allclose(data_point["target_forces"], [[0.1, 0.0, 0.0]] * 4)


def test_training_engine_load_data_with_delta(mock_db):
    """
    Tests the data loading logic when Delta Learning is enabled.
    """
    config = TrainingConfig(
        model_type="mace", learning_rate=0.01, num_epochs=1,
        r_cut=5.0, delta_learn=True, baseline_potential="lj"
    )
    engine = TrainingEngine(config, mock_db)

    # Get the expected baseline values
    row = mock_db.get(1)
    atoms = row.toatoms()
    expected_lj_energy, expected_lj_forces = calculate_lj_potential(atoms)

    prepared_data = engine._load_and_prepare_data(ids=[1])

    assert len(prepared_data) == 1
    data_point = prepared_data[0]

    # With delta learning, targets should be the DFT values minus the baseline
    expected_delta_energy = -1.0 - expected_lj_energy
    expected_delta_forces = np.array([[0.1, 0.0, 0.0]] * 4) - expected_lj_forces

    assert data_point["target_energy"] == pytest.approx(expected_delta_energy)
    assert np.allclose(data_point["target_forces"], expected_delta_forces)


def test_training_engine_skips_failed_calculations(mock_db):
    """
    Tests that the engine correctly skips database entries from failed calculations.
    """
    # Modify the mock to return a failed calculation
    mock_row = mock_db.get(1)
    mock_row.key_value_pairs["was_successful"] = False

    config = TrainingConfig(
        model_type="mace", learning_rate=0.01, num_epochs=1,
        r_cut=5.0, delta_learn=True, baseline_potential="lj"
    )
    engine = TrainingEngine(config, mock_db)

    prepared_data = engine._load_and_prepare_data(ids=[1])

    assert len(prepared_data) == 0
