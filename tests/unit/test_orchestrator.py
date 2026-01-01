from unittest.mock import MagicMock

import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTResults
from mlip_autopipec.orchestrator import Orchestrator


@pytest.fixture
def mock_db_wrapper() -> MagicMock:
    """Provides a MagicMock for the AseDBWrapper."""
    return MagicMock()


@pytest.fixture
def mock_structure_generator() -> MagicMock:
    """Provides a MagicMock for the StructureGenerator."""
    return MagicMock()


@pytest.fixture
def mock_labeling_engine() -> MagicMock:
    """Provides a MagicMock for the LabelingEngine."""
    return MagicMock()


@pytest.fixture
def mock_training_engine() -> MagicMock:
    """Provides a MagicMock for the TrainingEngine."""
    return MagicMock()


@pytest.fixture
def orchestrator(
    mock_db_wrapper: MagicMock,
    mock_structure_generator: MagicMock,
    mock_labeling_engine: MagicMock,
    mock_training_engine: MagicMock,
) -> Orchestrator:
    """Provides an Orchestrator instance with mocked dependencies."""
    return Orchestrator(
        db_wrapper=mock_db_wrapper,
        structure_generator=mock_structure_generator,
        labeling_engine=mock_labeling_engine,
        training_engine=mock_training_engine,
    )


def test_orchestrator_run_full_pipeline_success(
    orchestrator: Orchestrator,
    mock_db_wrapper: MagicMock,
    mock_structure_generator: MagicMock,
    mock_labeling_engine: MagicMock,
    mock_training_engine: MagicMock,
):
    """
    Tests the full 'Generate -> Label -> Train' workflow under ideal conditions.
    """
    # Arrange
    atoms_to_label = Atoms("Si")
    row_to_label = MagicMock(id=1, toatoms=MagicMock(return_value=atoms_to_label))
    mock_db_wrapper.get_rows_to_label.return_value = [row_to_label]

    dft_results = DFTResults(energy=-1.0, forces=[[0, 0, 0]], stress=[[0] * 3] * 3)
    mock_labeling_engine.execute.return_value = [(1, dft_results)]

    labeled_atoms = Atoms("Si", info={"energy": -1.0})
    labeled_row = MagicMock(toatoms=MagicMock(return_value=labeled_atoms))
    mock_db_wrapper.get_all_labeled_rows.return_value = [labeled_row]

    # Act
    orchestrator.run_full_pipeline()

    # Assert
    # Generation, Labeling, and Training should all be called in order.
    mock_structure_generator.execute.assert_called_once()
    mock_db_wrapper.get_rows_to_label.assert_called_once()
    mock_labeling_engine.execute.assert_called_once_with([(1, atoms_to_label)])
    mock_db_wrapper.update_row_with_dft_results.assert_called_once_with(1, dft_results)
    mock_db_wrapper.get_all_labeled_rows.assert_called_once()
    mock_training_engine.execute.assert_called_once_with([labeled_atoms])


def test_orchestrator_skips_labeling_if_no_new_structures(
    orchestrator: Orchestrator,
    mock_db_wrapper: MagicMock,
    mock_labeling_engine: MagicMock,
    mock_training_engine: MagicMock,
):
    """
    Tests that the labeling step is skipped if the DB has no unlabeled rows.
    """
    # Arrange
    mock_db_wrapper.get_rows_to_label.return_value = []  # No structures to label
    mock_db_wrapper.get_all_labeled_rows.return_value = [MagicMock()]  # But old ones to train

    # Act
    orchestrator.run_full_pipeline()

    # Assert
    mock_labeling_engine.execute.assert_not_called()
    mock_training_engine.execute.assert_called_once()  # Training should still run
