from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

# Under test
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.configs.models import MLIPTrainingConfig


@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper with labeled data."""
    db = MagicMock()

    # Create mock atoms with pre-calculated data
    atoms1 = Atoms("H", positions=[[0, 0, 0]])
    atoms1.calc = SinglePointCalculator(atoms1, energy=-1.0, forces=[[0.1, 0, 0]])

    atoms2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    atoms2.calc = SinglePointCalculator(atoms2, energy=-2.5, forces=[[0, 0, 0.2], [0, 0, -0.2]])

    mock_row1 = MagicMock()
    mock_row1.toatoms.return_value = atoms1
    mock_row2 = MagicMock()
    mock_row2.toatoms.return_value = atoms2

    db.get_labeled_rows.return_value = [mock_row1, mock_row2]
    return db


@pytest.fixture
def training_config():
    """Fixture for a simple MLIPTrainingConfig."""
    return MLIPTrainingConfig(
        model_type="ace",
        r_cut=4.5,
        delta_learning=False, # Test without delta learning first
        base_potential="lj_auto",
        loss_weights={"energy": 1.0, "force": 100.0},
    )


# Patch the actual training call, which is an external dependency
@patch("mlip_autopipec.modules.d_training_engine.train_ace_model")
def test_training_engine_run_no_delta(
    mock_train_call, mock_db_wrapper, training_config
):
    """
    Test the TrainingEngine's main run method without delta learning.
    """
    # --- Instantiate and Run ---
    engine = TrainingEngine(config=training_config, db_wrapper=mock_db_wrapper)
    engine.run()

    # --- Assertions ---
    # 1. Verify it fetched labeled data from the database
    mock_db_wrapper.get_labeled_rows.assert_called_once()

    # 2. Verify the training function was called
    mock_train_call.assert_called_once()

    # 3. Inspect the data passed to the training function
    call_args = mock_train_call.call_args[0]
    training_data = call_args[0]

    # Check that the data is a list of Atoms objects
    assert isinstance(training_data, list)
    assert len(training_data) == 2
    assert isinstance(training_data[0], Atoms)

    # Check that the energy and forces are the original DFT values
    assert training_data[0].get_potential_energy() == -1.0
    assert training_data[1].get_forces()[0][2] == 0.2


@patch("mlip_autopipec.modules.d_training_engine.train_ace_model")
@patch("mlip_autopipec.modules.d_training_engine.get_lj_potential")
def test_training_engine_run_with_delta(
    mock_get_lj, mock_train_call, mock_db_wrapper, training_config
):
    """
    Test the TrainingEngine's main run method with delta learning enabled.
    """
    # --- Setup Mocks ---
    training_config.delta_learning = True

    # Mock the Lennard-Jones potential to return predictable energy/forces
    mock_lj_calc = MagicMock()
    # Let's say LJ energy is -0.5 for H and -1.0 for H2
    # And forces are simple vectors
    mock_lj_calc.get_potential_energy.side_effect = [-0.5, -1.0]
    mock_lj_calc.get_forces.side_effect = [[[0.05, 0, 0]], [[0, 0, 0.05], [0, 0, -0.05]]]
    mock_get_lj.return_value = mock_lj_calc

    # --- Instantiate and Run ---
    engine = TrainingEngine(config=training_config, db_wrapper=mock_db_wrapper)
    engine.run()

    # --- Assertions ---
    mock_db_wrapper.get_labeled_rows.assert_called_once()
    mock_train_call.assert_called_once()

    # 3. Inspect the data passed to the training function - THIS IS THE KEY TEST
    call_args = mock_train_call.call_args[0]
    training_data = call_args[0]

    # Check that the energy and forces are now the DELTA values
    # DFT_energy - LJ_energy = -1.0 - (-0.5) = -0.5
    assert training_data[0].get_potential_energy() == pytest.approx(-0.5)

    # DFT_forces - LJ_forces = [0.1, 0, 0] - [0.05, 0, 0] = [0.05, 0, 0]
    assert training_data[0].get_forces()[0][0] == pytest.approx(0.05)

    # DFT_energy - LJ_energy = -2.5 - (-1.0) = -1.5
    assert training_data[1].get_potential_energy() == pytest.approx(-1.5)

    # DFT_forces - LJ_forces = [0, 0, 0.2] - [0, 0, 0.05] = [0, 0, 0.15]
    assert training_data[1].get_forces()[0][2] == pytest.approx(0.15)
