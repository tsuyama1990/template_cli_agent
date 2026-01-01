from unittest.mock import MagicMock

import pytest

from mlip_autopipec.config import MLIPTrainingConfig, ModelType
from mlip_autopipec.modules.training_engine import TrainingEngine


@pytest.fixture
def training_config():
    """Fixture for a sample MLIPTrainingConfig."""
    return MLIPTrainingConfig(
        model_type=ModelType.ACE,
        r_cut=5.0,
        loss_weights={"energy": 1.0, "forces": 100.0},
    )


def test_train_no_data(training_config):
    """Tests that the training engine handles the case where there is no data."""
    # Arrange
    engine = TrainingEngine(mlip_training_configuration=training_config)

    # Act
    engine.train([])

    # Assert
    # No assertion needed, just checking that it doesn't raise an exception
    pass


def test_train_raises_not_implemented_error(training_config):
    """Tests that the train method raises NotImplementedError."""
    # Arrange
    engine = TrainingEngine(mlip_training_configuration=training_config)

    # Act & Assert
    with pytest.raises(NotImplementedError):
        engine.train([(MagicMock(), MagicMock())])
