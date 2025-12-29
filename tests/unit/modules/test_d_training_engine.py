import os
from unittest.mock import Mock

import pytest

from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine

# --- Fixtures ---

@pytest.fixture
def mock_db_with_successful_data():
    """Mocks AseDB with a successful entry."""
    db = Mock()
    db.get.return_value = (Mock(), {"was_successful": True})
    return db

@pytest.fixture
def mock_db_with_failed_data():
    """Mocks AseDB with only a failed entry."""
    db = Mock()
    db.get.return_value = (Mock(), {"was_successful": False})
    return db

@pytest.fixture
def training_config():
    """Provides a standard TrainingConfig for testing."""
    return TrainingConfig(
        model_type="MACE",
        learning_rate=0.01,
        num_epochs=2,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="LJ"
    )

# --- Tests ---

def test_placeholder_engine_creates_file(training_config, mock_db_with_successful_data, tmp_path):
    """
    Tests that the placeholder TrainingEngine successfully creates a model file.
    """
    # 1. Setup
    output_dir = tmp_path / "models"

    # 2. Execute
    engine = TrainingEngine(config=training_config, db=mock_db_with_successful_data)
    model_path = engine.execute(ids=[1], output_dir=str(output_dir))

    # 3. Verify
    assert "mace_model_placeholder.pt" in model_path
    assert os.path.exists(model_path)
    with open(model_path) as f:
        content = f.read()
        assert "This is a placeholder model" in content

def test_placeholder_engine_raises_error_for_no_valid_data(training_config, mock_db_with_failed_data):
    """
    Tests that the engine raises a ValueError if there are no successful
    DFT calculations to "train" on.
    """
    # 1. Setup
    engine = TrainingEngine(config=training_config, db=mock_db_with_failed_data)

    # 2. Execute & Verify
    with pytest.raises(ValueError, match="No successful DFT calculations found"):
        engine.execute(ids=[1])
