# tests/unit/test_training_engine.py

from unittest.mock import MagicMock

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.configs.models import MLIPTrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db_wrapper_labeled():
    """Fixture to create a mock AseDBWrapper with labeled data."""
    db_wrapper = MagicMock()
    atoms = Atoms("H", cell=[1, 1, 1], pbc=True)
    kvp = {
        "id": 1,
        "state": "labeled",
        "data": {
            "energy": -10.0,
            "forces": [[0.1, 0.2, 0.3]],
        },
    }
    db_wrapper.get_labeled_atoms.return_value = [(atoms, kvp)]
    return db_wrapper


@pytest.fixture
def training_config_delta():
    """Fixture for a training config with delta learning enabled."""
    return MLIPTrainingConfig(
        model_type="ace",
        r_cut=5.0,
        delta_learning=True,
        base_potential="lj_auto",
        loss_weights={"energy": 1.0, "force": 100.0},
    )


@pytest.fixture
def training_config_no_delta():
    """Fixture for a training config with delta learning disabled."""
    return MLIPTrainingConfig(
        model_type="ace",
        r_cut=5.0,
        delta_learning=False,
        base_potential="lj_auto",
        loss_weights={"energy": 1.0, "force": 100.0},
    )


def test_training_engine_run_with_delta_learning(
    mocker, mock_db_wrapper_labeled, training_config_delta
):
    """Tests that the training engine correctly applies delta learning."""
    # Mock the baseline calculation to return predictable values
    mocker.patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._calculate_baseline",
        return_value={"energy": -2.0, "forces": np.array([[0.05, 0.05, 0.05]])},
    )
    # Mock the actual training call to inspect the data passed to it
    mock_training_call = mocker.patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._run_actual_training"
    )

    engine = TrainingEngine(
        config=training_config_delta, db_wrapper=mock_db_wrapper_labeled
    )
    engine.run()

    mock_training_call.assert_called_once()
    training_data = mock_training_call.call_args[0][0]
    assert len(training_data) == 1

    # Verify that the delta calculation is correct
    processed_atoms = training_data[0]
    assert processed_atoms.info["target_energy"] == pytest.approx(-10.0 - (-2.0))
    expected_forces = np.array([[0.1, 0.2, 0.3]]) - np.array([[0.05, 0.05, 0.05]])
    assert np.allclose(processed_atoms.info["target_forces"], expected_forces)


def test_training_engine_run_without_delta_learning(
    mocker, mock_db_wrapper_labeled, training_config_no_delta
):
    """Tests that the training engine runs correctly without delta learning."""
    mock_baseline = mocker.patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._calculate_baseline"
    )
    mock_training_call = mocker.patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._run_actual_training"
    )

    engine = TrainingEngine(
        config=training_config_no_delta, db_wrapper=mock_db_wrapper_labeled
    )
    engine.run()

    mock_baseline.assert_not_called()
    mock_training_call.assert_called_once()
    training_data = mock_training_call.call_args[0][0]
    assert len(training_data) == 1

    # Verify the target values are the original DFT values
    processed_atoms = training_data[0]
    assert processed_atoms.info["target_energy"] == -10.0
    assert np.allclose(
        processed_atoms.info["target_forces"], np.array([[0.1, 0.2, 0.3]])
    )


def test_training_engine_handles_no_data(mocker, training_config_delta):
    """Tests that the engine exits gracefully if no labeled data is found."""
    db_wrapper = MagicMock()
    db_wrapper.get_labeled_atoms.return_value = []
    mock_training_call = mocker.patch(
        "mlip_autopipec.modules.d_training_engine.TrainingEngine._run_actual_training"
    )

    engine = TrainingEngine(config=training_config_delta, db_wrapper=db_wrapper)
    engine.run()

    mock_training_call.assert_not_called()
