from unittest.mock import patch

import pytest
import yaml
from click.testing import CliRunner

from mlip_autopipec.cli import cli
from mlip_autopipec.data.models import Cycle01Config


@pytest.fixture
def runner():
    return CliRunner()

@pytest.fixture
def sample_config_file(tmp_path):
    """Creates a temporary valid config YAML file."""
    config_content = {
        'database_path': 'test.db',
        'dft_compute': {
            'code': 'quantum_espresso',
            'command': 'pw.x',
            'pseudopotentials': {
                'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
            },
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
    config_file = tmp_path / "config.yaml"
    with open(config_file, 'w') as f:
        yaml.dump(config_content, f)
    # Also create the dummy database file so FilePath validation passes
    (tmp_path / "test.db").touch()
    return str(config_file)


@patch('mlip_autopipec.workflow.Orchestrator')
def test_cli_run_cycle_success(mock_orchestrator, runner, sample_config_file):
    """Test the CLI `run-cycle` command with a valid config."""
    # Arrange
    mock_orchestrator_instance = mock_orchestrator.return_value
    mock_orchestrator_instance.run_label_and_train_workflow.return_value = None

    # Act
    result = runner.invoke(cli, ['run-cycle', '--config', sample_config_file])

    # Assert
    assert result.exit_code == 0

    # Check that Orchestrator was initialized with the correct config
    mock_orchestrator.assert_called_once()
    call_args = mock_orchestrator.call_args[0]
    assert len(call_args) == 1
    assert isinstance(call_args[0], Cycle01Config)

    # Check that the main workflow method was called
    mock_orchestrator_instance.run_label_and_train_workflow.assert_called_once()

def test_cli_run_cycle_file_not_found(runner):
    """Test the CLI with a non-existent config file."""
    result = runner.invoke(cli, ['run-cycle', '--config', 'nonexistent.yaml'])
    assert result.exit_code != 0
    assert "Path 'nonexistent.yaml' does not exist" in result.output

def test_cli_run_cycle_invalid_config(runner, tmp_path):
    """Test the CLI with an invalid config file."""
    invalid_config_content = {'database_path': 'test.db'} # Missing required fields
    config_file = tmp_path / "invalid_config.yaml"
    with open(config_file, 'w') as f:
        yaml.dump(invalid_config_content, f)
    (tmp_path / "test.db").touch()

    result = runner.invoke(cli, ['run-cycle', '--config', str(config_file)])
    assert result.exit_code != 0
    assert "Error loading configuration" in result.output
