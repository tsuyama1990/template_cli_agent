"""Unit tests for the CLI's core logic."""
from unittest.mock import patch, MagicMock
import pytest
from pydantic import ValidationError
import typer

# Import the functions and app object to be tested
from mlip_autopipec.cli import load_config, run
from mlip_autopipec.core.orchestrator import PipelineRunner
from mlip_autopipec.config.models import FullConfig


VALID_CONFIG_YAML = """
project_name: cli_test
system:
  elements: ['A', 'B']
  composition: {'A': 0.5, 'B': 0.5}
  lattice: 'fcc'
  num_structures: 1
exploration:
  temperature: 100.0
sampling:
  method: 'random'
  fraction: 1.0
"""

INVALID_CONFIG_YAML = """
project_name: cli_test_invalid
system:
  elements: ['A', 'B']
  composition: {'A': 0.6, 'B': 0.5} # Sums to 1.1
  lattice: 'fcc'
  num_structures: 1
exploration:
  temperature: 100.0
sampling:
  method: 'random'
  fraction: 1.0
"""

# --- Tests for the load_config function ---

def test_load_config_success(tmp_path):
    """Verify that load_config correctly parses a valid YAML file."""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(VALID_CONFIG_YAML)

    config = load_config(str(config_file))
    assert isinstance(config, FullConfig)
    assert config.project_name == "cli_test"

def test_load_config_file_not_found():
    """Verify that load_config raises typer.Exit when the file is not found."""
    with pytest.raises(typer.Exit) as e:
        load_config("non_existent_file.yaml")
    assert e.value.exit_code == 1

def test_load_config_validation_error(tmp_path):
    """Verify that load_config raises typer.Exit on a validation error."""
    config_file = tmp_path / "invalid_config.yaml"
    config_file.write_text(INVALID_CONFIG_YAML)

    with pytest.raises(typer.Exit) as e:
        load_config(str(config_file))
    assert e.value.exit_code == 1

# --- Test for the `run` command's logic by calling the function directly ---

@patch('mlip_autopipec.cli.load_config')
@patch('mlip_autopipec.cli.PipelineRunner')
def test_cli_run_function_logic(mock_pipeline_runner: MagicMock, mock_load_config: MagicMock):
    """
    Test the logic of the `run` function directly, bypassing the CliRunner.
    """
    # Arrange: Create a mock config object and a dummy config path
    mock_config = MagicMock(spec=FullConfig)
    mock_load_config.return_value = mock_config
    dummy_config_path = "dummy_path.yaml"

    # Act: Call the `run` function directly
    run(config=dummy_config_path)

    # Assert
    mock_load_config.assert_called_once_with(dummy_config_path)
    mock_pipeline_runner.assert_called_once_with(config=mock_config)
    mock_pipeline_runner.return_value.run.assert_called_once()
