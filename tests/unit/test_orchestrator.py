from unittest.mock import MagicMock

import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTResults
from mlip_autopipec.orchestrator import Orchestrator


@pytest.fixture
def mock_db_wrapper():
    """Provides a MagicMock for the AseDBWrapper."""
    return MagicMock()


@pytest.fixture
def mock_labeling_engine():
    """Provides a MagicMock for the LabelingEngine."""
    return MagicMock()


@pytest.fixture
def mock_training_engine():
    """Provides a MagicMock for the TrainingEngine."""
    return MagicMock()


@pytest.fixture
def orchestrator(
    mock_db_wrapper, mock_labeling_engine, mock_training_engine
) -> Orchestrator:
    """Provides an Orchestrator instance with mocked dependencies."""
    return Orchestrator(
        db_wrapper=mock_db_wrapper,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
    )


def test_orchestrator_run_workflow_success(
    orchestrator, mock_db_wrapper, mock_labeling_engine, mock_training_engine
):
    """
    Tests the full workflow execution of the Orchestrator under ideal conditions.
    """
    # Arrange
    # Mock data for labeling
    atoms_to_label = Atoms('Si')
    row_to_label = MagicMock()
    row_to_label.id = 1
    row_to_label.toatoms.return_value = atoms_to_label
    mock_db_wrapper.get_rows_to_label.return_value = [row_to_label]

    dft_results = DFTResults(energy=-1.0, forces=[[0, 0, 0]], stress=[[0]*3]*3)
    mock_labeling_engine.execute.return_value = [(1, dft_results)]

    # Mock data for training
    labeled_atoms = Atoms('Si', info={'energy': -1.0})
    labeled_row = MagicMock()
    labeled_row.toatoms.return_value = labeled_atoms
    mock_db_wrapper.get_all_labeled_rows.return_value = [labeled_row]

    # Act
    orchestrator.run_label_and_train_workflow()

    # Assert
    # Labeling assertions
    mock_db_wrapper.get_rows_to_label.assert_called_once()
    mock_labeling_engine.execute.assert_called_once_with([(1, atoms_to_label)])
    mock_db_wrapper.update_row_with_dft_results.assert_called_once_with(1, dft_results)

    # Training assertions
    mock_db_wrapper.get_all_labeled_rows.assert_called_once()
    mock_training_engine.execute.assert_called_once_with([labeled_atoms])


def test_orchestrator_no_rows_to_label(
    orchestrator, mock_db_wrapper, mock_labeling_engine, mock_training_engine
):
    """
    Tests that the labeling step is skipped when the database has no structures to label.
    """
    # Arrange
    mock_db_wrapper.get_rows_to_label.return_value = []
    mock_db_wrapper.get_all_labeled_rows.return_value = [MagicMock()] # Still train

    # Act
    orchestrator.run_label_and_train_workflow()

    # Assert
    mock_db_wrapper.get_rows_to_label.assert_called_once()
    mock_labeling_engine.execute.assert_not_called()
    mock_db_wrapper.update_row_with_dft_results.assert_not_called()
    mock_training_engine.execute.assert_called_once() # Training should still run


def test_orchestrator_no_labeled_rows_for_training(
    orchestrator, mock_db_wrapper, mock_labeling_engine, mock_training_engine
):
    """
    Tests that the training step is skipped when there are no labeled structures.
    """
    # Arrange
    mock_db_wrapper.get_rows_to_label.return_value = [MagicMock()] # Label first
    mock_labeling_engine.execute.return_value = [] # But labeling returns nothing
    mock_db_wrapper.get_all_labeled_rows.return_value = []

    # Act
    orchestrator.run_label_and_train_workflow()

    # Assert
    mock_labeling_engine.execute.assert_called_once()
    mock_db_wrapper.get_all_labeled_rows.assert_called_once()
    mock_training_engine.execute.assert_not_called()


def test_orchestrator_labeling_engine_raises_exception(
    orchestrator, mock_db_wrapper, mock_labeling_engine
):
    """
    Tests that the workflow halts and raises an exception if the labeling engine fails.
    """
    # Arrange
    mock_db_wrapper.get_rows_to_label.return_value = [MagicMock()]
    mock_labeling_engine.execute.side_effect = ValueError("DFT calculation failed")

    # Act & Assert
    with pytest.raises(ValueError, match="DFT calculation failed"):
        orchestrator.run_label_and_train_workflow()

    # Ensure training was not attempted
    mock_db_wrapper.get_all_labeled_rows.assert_not_called()
