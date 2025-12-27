from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.utils.baseline_potentials import calculate_lj_potential


@pytest.fixture
def h2o_atoms():
    """Returns a water molecule Atoms object."""
    from ase import Atoms

    return Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])


@pytest.fixture
def mock_db_with_data(h2o_atoms):
    """Creates a mock AseDB with a sample data row."""
    mock_db = MagicMock(spec=AseDB)
    mock_row = MagicMock()
    mock_row.toatoms.return_value = h2o_atoms
    mock_row.key_value_pairs = {
        "total_energy_ev": -450.0,
        "forces": np.random.rand(3, 3).tolist(),
        "stress": np.zeros((3, 3)).tolist(),
        "was_successful": True,
    }

    # Mock the connection attribute and its get method
    mock_db.connection.get.return_value = mock_row

    # Mock the _get_z_table method's dependency on the db
    with patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._get_z_table",
        return_value=MagicMock(),
    ):
        yield mock_db


def test_calculate_lj_potential(h2o_atoms):
    """Test the LJ potential calculation."""
    energy, forces = calculate_lj_potential(h2o_atoms)

    assert isinstance(energy, float)
    assert isinstance(forces, np.ndarray)
    assert forces.shape == (3, 3)
    # Check that forces sum to zero (Newton's third law)
    assert np.allclose(np.sum(forces, axis=0), 0)


@patch("mlip_autopipec.modules.d_training_engine.config_from_atoms")
def test_training_engine_load_data_delta_learn(mock_config_from_atoms, mock_db_with_data):
    """Test the TrainingEngine's data loading with Delta Learning."""
    config = TrainingConfig(
        model_type="mace",
        learning_rate=1e-3,
        num_epochs=10,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lj",
    )
    engine = TrainingEngine(config=config, db=mock_db_with_data)
    engine._load_and_prepare_data(ids=[1])

    mock_config_from_atoms.assert_called_once()
    atoms_arg = mock_config_from_atoms.call_args[0][0]

    assert atoms_arg.calc is not None
    dft_energy = mock_db_with_data.connection.get.return_value.key_value_pairs["total_energy_ev"]

    assert atoms_arg.get_potential_energy() != dft_energy
    assert atoms_arg.get_potential_energy() < dft_energy


@patch("mlip_autopipec.modules.d_training_engine.config_from_atoms")
def test_training_engine_load_data_no_delta(mock_config_from_atoms, mock_db_with_data):
    """Test the TrainingEngine's data loading without Delta Learning."""
    config = TrainingConfig(
        model_type="mace",
        learning_rate=1e-3,
        num_epochs=10,
        r_cut=5.0,
        delta_learn=False,
        baseline_potential="lj",
    )
    engine = TrainingEngine(config=config, db=mock_db_with_data)
    engine._load_and_prepare_data(ids=[1])

    mock_config_from_atoms.assert_called_once()
    atoms_arg = mock_config_from_atoms.call_args[0][0]

    assert atoms_arg.calc is not None
    dft_energy = mock_db_with_data.connection.get.return_value.key_value_pairs["total_energy_ev"]

    assert abs(atoms_arg.get_potential_energy() - dft_energy) < 1e-9
