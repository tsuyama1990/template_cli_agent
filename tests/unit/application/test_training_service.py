"""Unit tests for the TrainingService."""

import pytest
from unittest.mock import MagicMock, patch, ANY
from ase import Atoms
import numpy as np
import torch

from mlip_autopipec.application.services import TrainingService
from mlip_autopipec.config import Settings

@pytest.fixture
def mock_db_port():
    db = MagicMock()
    ar_dimer = Atoms('Ar2', positions=[[0, 0, 0], [0, 0, 3.8]], pbc=True, cell=np.diag([10, 10, 10]))
    ar_dimer.get_potential_energy = MagicMock(return_value=-0.02)
    ar_dimer.get_forces = MagicMock(return_value=np.array([[0, 0, 0.01], [0, 0, -0.01]]))
    db.get_atoms_by_state.return_value = [ar_dimer]
    return db

@pytest.fixture
def mock_settings():
    return Settings(model_path="test.pt", cutoff=5.0)

@patch('mlip_autopipec.application.services.torch.save')
@patch('mlip_autopipec.application.services.MACE')
@patch('mlip_autopipec.application.services.get_data_loader')
@patch('mlip_autopipec.application.services.WeightedEnergyForcesLoss')
@patch('mlip_autopipec.application.services.torch.optim.Adam')
def test_training_service_run(mock_adam, mock_loss, mock_loader, mock_mace, mock_torch_save, mock_db_port, mock_settings):
    # Arrange
    # Mock the data loader to return a single batch
    mock_loader.return_value = [{"some_data": "data"}]
    # Mock the loss function to return a mock loss object
    mock_loss_instance = MagicMock()
    mock_loss.return_value = mock_loss_instance
    mock_loss_value = MagicMock()
    mock_loss_instance.return_value = mock_loss_value

    # Mock the model to return a mock output
    mock_model_instance = MagicMock()
    mock_mace.return_value = mock_model_instance
    mock_model_instance.return_value = {"energy": torch.tensor(0.0)}

    service = TrainingService(mock_db_port, mock_settings)

    # Act
    service.run()

    # Assert
    mock_db_port.get_atoms_by_state.assert_called_once_with('labelled')
    mock_mace.assert_called_once()
    mock_loader.assert_called_once()
    mock_loss.assert_called_once()
    mock_adam.assert_called_once()

    # Verify that the core training steps were executed
    mock_loss_value.backward.assert_called()
    mock_adam.return_value.step.assert_called()

    mock_torch_save.assert_called_once()
