from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.data.models import Cycle01Config, DFTResults
from mlip_autopipec.orchestrator import Orchestrator


@pytest.fixture
def mock_cycle01_config(tmp_path):
    """Fixture to create a mock Cycle01Config object."""
    db_path = tmp_path / "test.db"
    db_path.touch()

    config_dict = {
        'database_path': db_path,
        'dft_compute': {
            'code': 'quantum_espresso',
            'command': 'pw.x',
            'pseudopotentials': {'Si': 'Si.UPF'},
            'ecutwfc': 60,
            'ecutrho': 240,
            'kpoints_density': 3,
        },
        'mlip_training': {
            'model_type': 'ace',
            'r_cut': 5.0,
            'delta_learning': False,
            'loss_weights': {'energy': 1.0, 'forces': 1.0}
        }
    }
    return Cycle01Config(**config_dict)


@patch('mlip_autopipec.orchestrator.TrainingEngine')
@patch('mlip_autopipec.orchestrator.LabelingEngine')
@patch('mlip_autopipec.orchestrator.AseDBWrapper')
def test_orchestrator_data_flow(
    mock_db_wrapper_cls, mock_labeling_engine_cls, mock_training_engine_cls,
    mock_cycle01_config
):
    """
    Tests that the Orchestrator correctly manages the data flow between the
    database and the engines.
    """
    # Arrange
    mock_db = MagicMock()
    mock_db_wrapper_cls.return_value = mock_db

    mock_row = MagicMock()
    mock_row.id = 1
    mock_row.toatoms.return_value = Atoms('Si')
    mock_db.get_rows_to_label.return_value = [mock_row]
    mock_db.get_all_labeled_rows.return_value = [mock_row]

    mock_labeler = MagicMock()
    mock_labeling_engine_cls.return_value = mock_labeler
    dft_results = DFTResults(energy=-1, forces=[[0,0,0]], stress=[0,0,0,0,0,0])
    mock_labeler.execute.return_value = [(1, dft_results)]

    mock_trainer = MagicMock()
    mock_training_engine_cls.return_value = mock_trainer

    # Act
    orchestrator = Orchestrator(mock_cycle01_config)
    orchestrator.run_label_and_train_workflow()

    # Assert
    mock_db.get_rows_to_label.assert_called_once()
    mock_labeler.execute.assert_called_once_with([(1, mock_row.toatoms())])
    mock_db.update_row_with_dft_results.assert_called_once_with(1, dft_results)
    mock_db.get_all_labeled_rows.assert_called_once()
    mock_trainer.execute.assert_called_once_with([mock_row.toatoms()])
