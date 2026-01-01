from unittest.mock import MagicMock

import pytest

from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import ILabelingEngine, IStructureGenerator, ITrainingEngine
from mlip_autopipec.workflow import WorkflowOrchestrator


@pytest.fixture
def mock_structure_generator():
    """Fixture for a mocked IStructureGenerator."""
    return MagicMock(spec=IStructureGenerator)


@pytest.fixture
def mock_labeling_engine():
    """Fixture for a mocked ILabelingEngine."""
    return MagicMock(spec=ILabelingEngine)


@pytest.fixture
def mock_training_engine():
    """Fixture for a mocked ITrainingEngine."""
    return MagicMock(spec=ITrainingEngine)


@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper."""
    return MagicMock(spec=AseDBWrapper)


def test_label_structure_by_id(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """Tests that label_structure_by_id correctly calls the labeling engine."""
    # Arrange
    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )
    structure_id = 123

    # Act
    orchestrator.label_structure_by_id(structure_id)

    # Assert
    mock_db_wrapper.get_atoms_by_id.assert_called_once_with(structure_id)
    mock_labeling_engine.label_structure.assert_called_once()
    mock_db_wrapper.update_labels.assert_called_once()


def test_label_structure_by_id_failure(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """Tests that label_structure_by_id correctly handles a labeling failure."""
    # Arrange
    mock_labeling_engine.label_structure.side_effect = Exception("Labeling failed")
    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )
    structure_id = 123

    # Act
    orchestrator.label_structure_by_id(structure_id)

    # Assert
    mock_db_wrapper.get_atoms_by_id.assert_called_once_with(structure_id)
    mock_labeling_engine.label_structure.assert_called_once()
    mock_db_wrapper.update_state.assert_called_once_with(structure_id, "labeling_failed")


def test_run_training(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """Tests that run_training correctly calls the training engine."""
    # Arrange
    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )

    # Act
    orchestrator.run_training()

    # Assert
    mock_db_wrapper.get_all_labeled_atoms.assert_called_once()
    mock_training_engine.train.assert_called_once()


def test_run_training_failure(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """Tests that run_training correctly handles a training failure."""
    # Arrange
    mock_training_engine.train.side_effect = Exception("Training failed")
    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )

    # Act
    orchestrator.run_training()

    # Assert
    mock_db_wrapper.get_all_labeled_atoms.assert_called_once()
    mock_training_engine.train.assert_called_once()
