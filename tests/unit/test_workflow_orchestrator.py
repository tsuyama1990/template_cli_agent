from unittest.mock import MagicMock

import pytest

from mlip_autopipec.interfaces import ILabelingEngine, ITrainingEngine
from mlip_autopipec.workflow import WorkflowOrchestrator


@pytest.fixture
def mock_labeling_engine():
    """Fixture for a mocked ILabelingEngine."""
    return MagicMock(spec=ILabelingEngine)


@pytest.fixture
def mock_training_engine():
    """Fixture for a mocked ITrainingEngine."""
    return MagicMock(spec=ITrainingEngine)


def test_label_structure_by_id(mock_labeling_engine, mock_training_engine):
    """Tests that label_structure_by_id correctly calls the labeling engine."""
    # Arrange
    orchestrator = WorkflowOrchestrator(
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
    )
    structure_id = 123

    # Act
    orchestrator.label_structure_by_id(structure_id)

    # Assert
    mock_labeling_engine.label_structure.assert_called_once_with(structure_id)


def test_run_training(mock_labeling_engine, mock_training_engine):
    """Tests that run_training correctly calls the training engine."""
    # Arrange
    orchestrator = WorkflowOrchestrator(
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
    )

    # Act
    orchestrator.run_training()

    # Assert
    mock_training_engine.train.assert_called_once()
