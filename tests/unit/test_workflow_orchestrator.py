from unittest.mock import MagicMock, call

import pytest
from ase import Atoms

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


def test_run_workflow_with_empty_db_generates_and_labels(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """
    Tests that the main run() method correctly calls the structure generator
    when the database is empty and then proceeds to label the new structures.
    """
    # Arrange
    mock_db_wrapper.is_empty.return_value = True
    generated_atoms = [Atoms("H"), Atoms("He")]
    mock_structure_generator.generate.return_value = generated_atoms
    mock_db_wrapper.get_unlabeled_ids.return_value = [1, 2] # IDs of the new atoms

    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )

    # Act
    orchestrator.run()

    # Assert
    mock_db_wrapper.is_empty.assert_called_once()
    mock_structure_generator.generate.assert_called_once()

    # Check that the orchestrator saved the atoms
    expected_calls = [
        call(generated_atoms[0], state="unlabeled"),
        call(generated_atoms[1], state="unlabeled"),
    ]
    mock_db_wrapper.add_atoms.assert_has_calls(expected_calls, any_order=True)

    # Check that it proceeded to label
    assert mock_db_wrapper.get_unlabeled_ids.called
    assert mock_labeling_engine.label_structure.call_count == 2

    # Check that it ran training at the end
    mock_training_engine.train.assert_called_once()


def test_run_workflow_with_populated_db_skips_generation(
    mock_structure_generator, mock_labeling_engine, mock_training_engine, mock_db_wrapper
):
    """
    Tests that the main run() method skips structure generation if the
    database is not empty.
    """
    # Arrange
    mock_db_wrapper.is_empty.return_value = False
    mock_db_wrapper.get_unlabeled_ids.return_value = [3, 4] # Pre-existing IDs

    orchestrator = WorkflowOrchestrator(
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
        db_wrapper=mock_db_wrapper,
    )

    # Act
    orchestrator.run()

    # Assert
    mock_db_wrapper.is_empty.assert_called_once()
    mock_structure_generator.generate.assert_not_called()
    mock_db_wrapper.add_atoms.assert_not_called()

    # Check that it proceeded to label the existing structures
    assert mock_labeling_engine.label_structure.call_count == 2
    mock_training_engine.train.assert_called_once()
