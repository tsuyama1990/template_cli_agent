from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.data.models import MLIPTraining
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db_wrapper_labeled():
    db = MagicMock()
    atoms1 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], cell=[10, 10, 10])
    atoms2 = Atoms("Si2", positions=[[0, 0, 0], [1, 1, 1]], cell=[10, 10, 10])

    row1 = MagicMock()
    row1.toatoms.return_value = atoms1
    row2 = MagicMock()
    row2.toatoms.return_value = atoms2

    db.connection.select.return_value = [row1, row2]
    return db


@pytest.fixture
def mock_training_config():
    return MLIPTraining(
        model_type="ace",
        r_cut=5.0,
        delta_learning=False,
        loss_weights={"energy": 1.0, "forces": 1.0},
    )


@patch("torch.save")
@patch("mlip_autopipec.modules.d_training_engine.MACE")
def test_training_engine_execute(
    mock_mace_model, mock_torch_save, mock_db_wrapper_labeled, mock_training_config
):
    # Arrange
    engine = TrainingEngine(mock_db_wrapper_labeled, mock_training_config)

    # Act
    engine.execute()

    # Assert
    # Check that data was queried
    mock_db_wrapper_labeled.connection.select.assert_called_once_with(labeled=True)

    # Check that model was built and saved
    mock_mace_model.assert_called_once()
    mock_torch_save.assert_called_once()

    # Check that the object passed to torch.save is an instance of the mocked MACE class
    saved_object = mock_torch_save.call_args[0][0]
    assert isinstance(saved_object, MagicMock)  # It's a mock of the MACE class instance


def test_training_engine_no_data(mock_db_wrapper_labeled, mock_training_config):
    # Arrange
    mock_db_wrapper_labeled.connection.select.return_value = []
    engine = TrainingEngine(mock_db_wrapper_labeled, mock_training_config)

    # Act
    with patch("builtins.print") as mock_print:
        engine.execute()

    # Assert
    mock_print.assert_called_with("No labeled structures found to train on.")
