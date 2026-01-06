"""Unit tests for the command-line interface."""

from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from mlip_autopipec.cli import app

runner = CliRunner()


@patch("mlip_autopipec.cli.PipelineRunner")
def test_cli_successful_invocation(mock_pipeline_runner, tmp_path: Path) -> None:
    # Corresponds to UAT-C1-001: Successful End-to-End Run with a Valid Configuration
    # This test simulates the primary "happy path" for the user.
    """Test a successful invocation of the CLI with a valid config."""
    # 1. Create a dummy valid config file
    config_content = """
    project_name: test_cli_success
    system:
      elements: ['Fe', 'Pt']
      composition: {'Fe': 0.5, 'Pt': 0.5}
      lattice: 'fcc'
      num_structures: 1
    exploration:
      temperature: 300.0
    sampling:
      method: 'random'
      fraction: 1.0
    """
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    # 2. Invoke the CLI
    result = runner.invoke(app, ["run", "--config", str(config_file)])

    # 3. Assertions
    assert result.exit_code == 0
    assert "MLIP-AutoPipe finished successfully!" in result.stdout
    mock_pipeline_runner.assert_called_once()
    mock_pipeline_runner.return_value.run.assert_called_once()


def test_cli_file_not_found() -> None:
    # Corresponds to UAT-C1-002: Application Handles a Missing Configuration File Gracefully
    # Verifies that the CLI provides a clear error for a common user mistake.
    """Test that the CLI handles a non-existent config file gracefully."""
    result = runner.invoke(app, ["run", "--config", "non_existent_file.yaml"])

    # Typer's `exists=True` handles this, so the error message is from Typer
    assert result.exit_code != 0
    assert "does not exist" in result.stdout


@patch("mlip_autopipec.cli.PipelineRunner")
def test_cli_invalid_config(mock_pipeline_runner, tmp_path: Path) -> None:
    # Corresponds to UAT-C1-003: Application Rejects a Configuration with Invalid Schema
    # Verifies that Pydantic validation errors are caught and reported to the user.
    """Test that the CLI handles a configuration with a validation error."""
    # 1. Create a config with a composition sum error
    config_content = """
    project_name: test_cli_invalid
    system:
      elements: ['Fe', 'Pt']
      composition: {'Fe': 0.6, 'Pt': 0.5} # Invalid
      lattice: 'fcc'
      num_structures: 1
    exploration:
      temperature: 300.0
    sampling:
      method: 'random'
      fraction: 1.0
    """
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    # 2. Invoke the CLI
    result = runner.invoke(app, ["run", "--config", str(config_file)])

    # 3. Assertions
    assert result.exit_code == 1
    assert "Configuration validation failed" in result.stdout
    mock_pipeline_runner.assert_not_called()  # Pipeline should not be instantiated
