from unittest.mock import MagicMock, call, patch

import pytest

from mlip_autopipec.data.models import Cycle01Config
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
def test_orchestrator_workflow_sequence(
    mock_labeling_engine, mock_training_engine, mock_cycle01_config
):
    """
    Tests that the Orchestrator initializes and runs the engines in the correct
    Label -> Train sequence.
    """
    mock_labeler_instance = MagicMock()
    mock_trainer_instance = MagicMock()
    mock_labeling_engine.return_value = mock_labeler_instance
    mock_training_engine.return_value = mock_trainer_instance

    manager = MagicMock()
    manager.attach_mock(mock_labeler_instance.execute, 'labeler_execute')
    manager.attach_mock(mock_trainer_instance.execute, 'trainer_execute')

    orchestrator = Orchestrator(mock_cycle01_config)
    orchestrator.run_label_and_train_workflow()

    mock_labeling_engine.assert_called_once()
    mock_training_engine.assert_called_once()

    expected_call_order = [
        call.labeler_execute(),
        call.trainer_execute()
    ]
    assert manager.mock_calls == expected_call_order
