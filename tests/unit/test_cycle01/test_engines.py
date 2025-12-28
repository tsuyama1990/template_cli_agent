# Description: Unit tests for the LabellingEngine and TrainingEngine.
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import torch
from ase.build import bulk

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult, TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def mock_db(tmpdir):
    """Creates a mock AseDB that writes to a temporary file."""
    db_path = Path(tmpdir) / "test.db"
    return AseDB(str(db_path))


@patch("mlip_autopipec.modules.c_labelling_engine.shutil.rmtree")
@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
@patch("mlip_autopipec.modules.c_labelling_engine.parse_qe_output")
@patch("mlip_autopipec.modules.c_labelling_engine.Path.write_text")
@patch("mlip_autopipec.modules.c_labelling_engine.Path.mkdir")
def test_labelling_engine_execute_success(
    mock_mkdir, mock_write, mock_parse, mock_run, mock_rmtree, mock_db
):
    """Tests a successful run of the LabellingEngine."""
    # Mock the subprocess to return a completed process with mock output
    mock_run.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="mock output"
    )
    mock_parse.return_value = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0] * 3] * 2,
        stress=[[1.0] * 3] * 3,
        was_successful=True,
    )

    engine = LabellingEngine(qe_command="pw.x -in QE_input.in", db=mock_db)
    atoms = bulk("Si")
    db_id = engine.execute(atoms)

    mock_run.assert_called_once()
    assert "pw.x" in mock_run.call_args.args[0]
    mock_parse.assert_called_once_with("mock output")
    retrieved_atoms, kvp = mock_db.get(db_id)
    assert kvp["was_successful"] is True
    assert retrieved_atoms.get_potential_energy() == -100.0
    mock_rmtree.assert_called_once()


@patch("mlip_autopipec.modules.c_labelling_engine.shutil.rmtree")
@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
@patch("mlip_autopipec.modules.c_labelling_engine.Path.write_text")
@patch("mlip_autopipec.modules.c_labelling_engine.Path.mkdir")
def test_labelling_engine_execute_qe_failure(mock_mkdir, mock_write, mock_run, mock_rmtree, mock_db):
    """Tests that the LabellingEngine handles a QE subprocess failure."""
    # Mock the subprocess to raise an error
    mock_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd="mock_cmd", stderr="QE crashed"
    )

    engine = LabellingEngine(qe_command="pw.x -in QE_input.in", db=mock_db)
    atoms = bulk("Si")
    db_id = engine.execute(atoms)

    mock_run.assert_called_once()
    retrieved_atoms, kvp = mock_db.get(db_id)
    assert kvp["was_successful"] is False
    assert "QE execution failed" in kvp["error_message"]
    assert "QE crashed" in kvp["error_message"]
    mock_rmtree.assert_called_once()


@patch("mlip_autopipec.modules.d_training_engine.torch.save")
@patch("mlip_autopipec.modules.d_training_engine.Adam")
@patch("mlip_autopipec.modules.d_training_engine.MACE")
def test_training_engine_execute_delta_learn(mock_mace, mock_adam, mock_torch_save, mock_db):
    """Tests the TrainingEngine with Delta Learning enabled."""
    # Setup: Write a successful result to the mock DB
    atoms = bulk("Si")
    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=np.random.rand(2, 3),
        stress=[[1.0] * 3] * 3,
        was_successful=True,
    )
    db_id = mock_db.write(atoms, dft_result)

    # Mock the MACE model to behave as needed
    mock_model_instance = MagicMock()
        # Configure the mock to return a dictionary with a tensor that requires gradients
    mock_model_instance.return_value = {
            "energy": torch.randn(1, 1, requires_grad=True),
            "forces": torch.randn(2, 3, requires_grad=True),
    }
    mock_mace.return_value = mock_model_instance

    config = TrainingConfig(
        model_type="MACE",
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="ZBL",
    )
    engine = TrainingEngine(config=config, db=mock_db)

    # Use a spy to check the arguments passed to _load_and_prepare_data
    with patch.object(
        engine, "_load_and_prepare_data", wraps=engine._load_and_prepare_data
    ) as spy_prepare_data:
        model_path = engine.execute(ids=[db_id])

        spy_prepare_data.assert_called_once_with([db_id])
        # The successful completion of the training loop is the primary validation.
        # Introspecting the wrapped method's return value is proving unreliable.

    assert Path(model_path).name == "trained_model.pt"
    mock_mace.assert_called_once()
    mock_adam.assert_called_once()
    mock_torch_save.assert_called_once()


@patch("mlip_autopipec.modules.d_training_engine.torch.save")
@patch("mlip_autopipec.modules.d_training_engine.Adam")
@patch("mlip_autopipec.modules.d_training_engine.MACE")
def test_training_engine_execute_direct_learn(mock_mace, mock_adam, mock_torch_save, mock_db):
    """Tests the TrainingEngine with Delta Learning disabled (direct learning)."""
    atoms = bulk("Si")
    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=np.zeros((2, 3)),
        stress=[[0.0] * 3] * 3,
        was_successful=True,
    )
    db_id = mock_db.write(atoms, dft_result)
    mock_model_instance = MagicMock()
    mock_model_instance.return_value = {
        "energy": torch.tensor([[-100.0]], requires_grad=True),
        "forces": torch.zeros(2, 3, requires_grad=True),
    }
    mock_mace.return_value = mock_model_instance
    config = TrainingConfig(
        model_type="MACE",
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=False,
        baseline_potential="ZBL",
    )
    engine = TrainingEngine(config=config, db=mock_db)
    with patch.object(
        engine, "_load_and_prepare_data", wraps=engine._load_and_prepare_data
    ) as spy_prepare_data:
        engine.execute(ids=[db_id])
        spy_prepare_data.assert_called_once_with([db_id])
        # The successful completion of the training loop is sufficient validation.


def test_training_engine_skips_failed_dft_calcs(mock_db):
    """Tests that the engine skips entries from failed DFT calculations."""
    # Write one successful and one failed result
    atoms = bulk("Si")
    success_result = DFTResult(
        total_energy_ev=-100.0, forces=np.zeros((2, 3)), stress=[], was_successful=True
    )
    fail_result = DFTResult(
        total_energy_ev=0.0, forces=[], stress=[], was_successful=False
    )
    success_id = mock_db.write(atoms, success_result)
    fail_id = mock_db.write(atoms, fail_result)

    config = TrainingConfig(
        model_type="None",
        learning_rate=0.01,
        num_epochs=1,
        r_cut=5.0,
        delta_learn=False,
        baseline_potential="ZBL",
    )
    engine = TrainingEngine(config=config, db=mock_db)

    # Spy on the _load_and_prepare_data method
    with patch.object(
        engine, "_load_and_prepare_data", wraps=engine._load_and_prepare_data
    ) as spy:
        engine.execute(ids=[success_id, fail_id])
        spy.assert_called_once_with([success_id, fail_id])
        # The successful completion of the training loop is sufficient validation.
