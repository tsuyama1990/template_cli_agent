# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import torch
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db_with_data() -> MagicMock:
    """Provides a mock AseDB instance with a sample data record."""
    db = MagicMock()
    atoms = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    # Attach a calculator with dummy results
    atoms.calc = SinglePointCalculator(atoms, energy=-10.0, forces=np.array([[0.1, 0.2, 0.3]]))
    db.get.return_value = (atoms, {"was_successful": True})
    return db


@pytest.fixture
def training_config() -> TrainingConfig:
    """Provides a sample TrainingConfig."""
    return TrainingConfig(
        model_type="mace",
        learning_rate=0.001,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lj",
    )


@patch("mlip_autopipec.modules.d_training_engine.MACE")
@patch("torch.save")
def test_training_engine_execute(
    mock_torch_save: MagicMock,
    mock_mace: MagicMock,
    mock_db_with_data: MagicMock,
    training_config: TrainingConfig,
):
    """Test the main execution flow of the TrainingEngine, mocking the ML model."""
    # Mock the MACE model to avoid actual training
    mock_model_instance = MagicMock()
    mock_model_instance.parameters.return_value = [
        torch.nn.Parameter(torch.randn(1, requires_grad=True))
    ]
    # The model's output must be a dict with tensors requiring gradients for loss.backward()
    mock_model_instance.return_value = {
        "energy": torch.tensor([0.0], requires_grad=True, dtype=torch.float64),
        "forces": torch.tensor([[0.0, 0.0, 0.0]], requires_grad=True, dtype=torch.float64),
    }
    mock_mace.return_value = mock_model_instance

    engine = TrainingEngine(config=training_config, db=mock_db_with_data)

    model_path = engine.execute(ids=[1])

    # --- Verifications ---
    # 1. Database was called
    mock_db_with_data.get.assert_called_once_with(1)

    # 2. MACE model was initialized
    mock_mace.assert_called_once()

    # 3. Training loop ran (model was called)
    assert mock_model_instance.call_count > 0

    # 4. torch.save was called with the model and a valid path
    mock_torch_save.assert_called_once()
    saved_model = mock_torch_save.call_args[0][0]
    saved_path = mock_torch_save.call_args[0][1]
    assert saved_model is mock_model_instance
    assert str(saved_path) == "models/trained_model.pt"

    # 5. A valid model path string was returned
    assert model_path == "models/trained_model.pt"
