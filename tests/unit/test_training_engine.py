import pytest
from unittest.mock import MagicMock, patch
import os
import json
import numpy as np

from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.config import MLIPTrainingConfig, DFTResult
from mlip_autopipec.database import AseDBWrapper
from ase.build import molecule

@pytest.fixture
def mock_db_wrapper():
    """Provides a mock AseDBWrapper."""
    return MagicMock(spec=AseDBWrapper)

@pytest.fixture
def training_config():
    """Provides a standard MLIPTrainingConfig."""
    return MLIPTrainingConfig(r_cut=5.0, loss_weights={"energy": 1.0})

def test_training_engine_no_data(mock_db_wrapper, training_config, capsys):
    """
    Tests that the training engine correctly handles the case where
    the database has no labeled data.
    """
    # Configure the mock to return an empty list
    mock_db_wrapper.get_all_labeled_atoms.return_value = []

    engine = TrainingEngine(mock_db_wrapper, training_config)
    engine.train()

    # Verify that the correct message was printed and no model was created
    captured = capsys.readouterr()
    assert "No labeled data found. Aborting training." in captured.out
    assert not os.path.exists("ace_model.json")

def test_training_engine_with_data(mock_db_wrapper, training_config, capsys):
    """
    Tests that the training engine successfully "trains" and creates a
    placeholder model file when labeled data is present.
    """
    # Create some mock labeled data
    h2o = molecule("H2O")
    dft_result = DFTResult(
        energy=-450.0,
        forces=np.zeros((3, 3)),
        stress=np.zeros((3, 3))
    )
    mock_labeled_data = [(h2o, dft_result)]

    # Configure the mock to return the data
    mock_db_wrapper.get_all_labeled_atoms.return_value = mock_labeled_data

    engine = TrainingEngine(mock_db_wrapper, training_config)
    engine.train()

    # Verify that the correct messages were printed
    captured = capsys.readouterr()
    assert "Found 1 labeled structures for training." in captured.out
    assert "Training complete. Model saved to 'ace_model.json'." in captured.out

    # Verify that the model file was created and contains the correct info
    model_filename = "ace_model.json"
    assert os.path.exists(model_filename)

    with open(model_filename, 'r') as f:
        model_data = json.load(f)
        assert model_data["training_set_size"] == 1
        assert model_data["config"]["r_cut"] == 5.0

    # Clean up the created file
    os.remove(model_filename)
