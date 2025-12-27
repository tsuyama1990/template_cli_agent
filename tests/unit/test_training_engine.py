import pytest
from unittest.mock import MagicMock, patch

from ase.build import bulk

from src.mlip_autopipec.data.database import AseDB
from src.mlip_autopipec.data.models import TrainingConfig, DFTResult
from src.mlip_autopipec.modules.d_training_engine import TrainingEngine

@pytest.fixture
def mock_db():
    db = MagicMock(spec=AseDB)
    si_atoms = bulk('Si', 'diamond', a=5.43)
    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]],
        stress=[[0.0]*3]*3,
        was_successful=True
    )
    failed_result = DFTResult(was_successful=False)

    def get_side_effect(db_id):
        if db_id == 1:
            return (si_atoms.copy(), dft_result)
        if db_id == 2:
            return (si_atoms.copy(), failed_result)
        return None

    db.get.side_effect = get_side_effect
    return db

@pytest.fixture
def training_config():
    return TrainingConfig(
        model_type="MACE", learning_rate=0.01, num_epochs=10,
        r_cut=5.0, delta_learn=True, baseline_potential="LJ"
    )

@patch('subprocess.run')
def test_execute_success(mock_subprocess_run, training_config, mock_db, tmp_path):
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_subprocess_run.return_value = mock_process

    engine = TrainingEngine(config=training_config, db=mock_db)
    model_save_path = tmp_path / "test_model.pt"

    result_path = engine.execute(ids=[1], model_save_path=str(model_save_path), work_dir=tmp_path / "mace_temp")

    assert result_path == str(model_save_path)
    mock_subprocess_run.assert_called_once()

    # Check that the command was constructed correctly
    args, kwargs = mock_subprocess_run.call_args
    command = args[0]
    assert "mace_run_train" in command
    assert f"--name=MLIP_AutoPipe_MACE" in command
    assert f"--train_file={tmp_path / 'mace_temp/train.xyz'}" in command
    assert f"--model_dir={tmp_path}" in command
    assert f"--num_epochs=10" in command

@patch('subprocess.run')
def test_execute_training_failure(mock_subprocess_run, training_config, mock_db, tmp_path):
    mock_process = MagicMock()
    mock_process.returncode = 1
    mock_process.stderr = "MACE training error"
    mock_subprocess_run.return_value = mock_process

    engine = TrainingEngine(config=training_config, db=mock_db)
    model_save_path = tmp_path / "test_model.pt"

    result_path = engine.execute(ids=[1], model_save_path=str(model_save_path), work_dir=tmp_path / "mace_temp")

    assert result_path is None

def test_execute_no_valid_data(training_config, mock_db, tmp_path):
    engine = TrainingEngine(config=training_config, db=mock_db)
    model_save_path = tmp_path / "test_model.pt"
    result_path = engine.execute(ids=[2, 3], model_save_path=str(model_save_path), work_dir=tmp_path / "mace_temp")

    assert result_path is None
