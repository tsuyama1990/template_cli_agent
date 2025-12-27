import pytest
from unittest.mock import MagicMock, patch
from click.testing import CliRunner
from pathlib import Path
import yaml

from ase.build import bulk

from src.mlip_autopipec.main_cycle01 import main
from src.mlip_autopipec.data.models import DFTResult, TrainingConfig

@pytest.fixture
def test_structure(tmp_path):
    """Creates a simple structure file for testing."""
    structure_file = tmp_path / "si.cif"
    atoms = bulk('Si', 'diamond', a=5.43)
    atoms.write(structure_file)
    return structure_file

@pytest.fixture
def test_config(tmp_path):
    """Creates a detailed config file for testing."""
    config_file = tmp_path / "config.yaml"
    config_data = {
        "db_path": str(tmp_path / "test.db"),
        "qe_command": "mpirun -np 4 pw.x",
        "training": {
            "model_type": "MACE",
            "learning_rate": 0.05,
            "num_epochs": 50,
            "r_cut": 4.5,
            "delta_learn": False,
            "baseline_potential": "LJ"
        }
    }
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f)
    return config_file, config_data

@patch('src.mlip_autopipec.orchestrator_cycle01.LabellingEngine')
@patch('src.mlip_autopipec.orchestrator_cycle01.TrainingEngine')
@patch('src.mlip_autopipec.orchestrator_cycle01.AseDB')
def test_e2e_successful_workflow(mock_asedb, mock_training_engine, mock_labelling_engine, test_config):
    config_file, config_data = test_config
    structure_file = Path(config_file).parent / "si.cif"
    bulk('Si', 'diamond', a=5.43).write(structure_file)

    mock_labeller_instance = MagicMock()
    mock_labeller_instance.execute.return_value = 1
    mock_labelling_engine.return_value = mock_labeller_instance

    mock_trainer_instance = MagicMock()
    mock_trainer_instance.execute.return_value = "trained_model.pt"
    mock_training_engine.return_value = mock_trainer_instance

    mock_db_instance = MagicMock()
    mock_db_instance.get.return_value = (None, DFTResult(was_successful=True))
    mock_asedb.return_value = mock_db_instance

    runner = CliRunner()
    result = runner.invoke(
        main,
        ['--config', str(config_file), '--structure', str(structure_file)]
    )

    assert result.exit_code == 0
    assert "Workflow complete" in result.output

    # Verify that the components were initialized with the correct config values
    mock_asedb.assert_called_once_with(Path(config_data["db_path"]))
    mock_labelling_engine.assert_called_once_with(
        qe_command=config_data["qe_command"],
        db=mock_db_instance
    )

    expected_training_config = TrainingConfig(**config_data["training"])
    mock_training_engine.assert_called_once_with(
        config=expected_training_config,
        db=mock_db_instance
    )
