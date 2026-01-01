from unittest.mock import MagicMock

from mlip_autopipec.data.models import MLIPTraining
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def test_training_engine_initialization():
    """Tests that the TrainingEngine can be initialized."""
    mock_config = MagicMock(spec=MLIPTraining)
    engine = TrainingEngine(config=mock_config)
    assert isinstance(engine, TrainingEngine)
