from unittest.mock import MagicMock

import pytest

from mlip_autopipec.config import MLIPTrainingConfig, ModelType
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.modules.training_engine import TrainingEngine


@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper."""
    return MagicMock(spec=AseDBWrapper)


@pytest.fixture
def training_config():
    """Fixture for a sample MLIPTrainingConfig."""
    return MLIPTrainingConfig(
        model_type=ModelType.ACE,
        r_cut=5.0,
        loss_weights={"energy": 1.0, "forces": 100.0},
    )


def test_train_no_data(mock_db_wrapper, training_config):
    """Tests that the training engine handles the case where there is no data."""
    # Arrange
    mock_db_wrapper.get_all_labeled_atoms.return_value = []
    engine = TrainingEngine(
        training_config=training_config, db_wrapper=mock_db_wrapper
    )

    # Act
    engine.train()

    # Assert
    mock_db_wrapper.get_all_labeled_atoms.assert_called_once()


def test_train_raises_not_implemented_error(mock_db_wrapper, training_config):
    """Tests that the train method raises NotImplementedError."""
    # Arrange
    mock_db_wrapper.get_all_labeled_atoms.return_value = [(MagicMock(), MagicMock())]
    engine = TrainingEngine(
        training_config=training_config, db_wrapper=mock_db_wrapper
    )

    # Act & Assert
    with pytest.raises(NotImplementedError):
        engine.train()
