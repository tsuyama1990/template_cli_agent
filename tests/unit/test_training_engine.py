from unittest.mock import MagicMock, patch

import pytest
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db():
    """Provides a mock AseDB instance."""
    db = MagicMock(spec=AseDB)
    atoms = bulk("Si", "diamond", a=5.43)
    db.get.return_value = (atoms, {"was_successful": True})
    return db


@pytest.fixture
def training_config():
    """Provides a default TrainingConfig."""
    return TrainingConfig(
        model_type="mace",
        learning_rate=0.01,
        num_epochs=10,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="lj",
    )


@patch("mlip_autopipec.modules.d_training_engine.ase_write")
@patch("mlip_autopipec.modules.d_training_engine.subprocess.run")
def test_training_engine_execute(mock_subprocess_run, mock_ase_write, mock_db, training_config):
    """
    Tests that the TrainingEngine correctly constructs the mace_run_train command
    and writes the training data.
    """
    engine = TrainingEngine(config=training_config, db=mock_db)
    model_path = "test_models/my_model.pt"

    engine.execute(ids=[1, 2], model_save_path=model_path)

    # 1. Verify that ase.io.write was called correctly
    mock_ase_write.assert_called_once()
    args, kwargs = mock_ase_write.call_args
    assert "train.xyz" in args[0]  # Check that the temp file is being written to
    assert len(args[1]) == 2  # Check that we are writing 2 atoms objects

    # 2. Verify that subprocess.run was called with the correct command
    mock_subprocess_run.assert_called_once()
    command = mock_subprocess_run.call_args[0][0]

    assert "mace_run_train" in command
    assert f"--name={model_path.split('/')[-1]}" in command
    assert f"--r_max={training_config.r_cut}" in command
    assert f"--max_num_epochs={training_config.num_epochs}" in command
    assert f"--results_dir={model_path.rsplit('/', 1)[0]}" in command
