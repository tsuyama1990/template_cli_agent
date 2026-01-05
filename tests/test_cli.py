"""Unit tests for the command-line interface."""

from pathlib import Path

from typer.testing import CliRunner

from mlip_autopipec.cli import app

runner = CliRunner()


def test_cli_successful_run(mocker, tmp_path: Path):
    """Test a successful run of the CLI."""
    mock_pipeline_runner = mocker.patch("mlip_autopipec.cli.PipelineRunner")

    config_content = """
    project_name: test_project
    system:
      elements: ['Fe', 'Pt']
      composition: {'Fe': 0.5, 'Pt': 0.5}
      lattice: 'fcc'
      num_structures: 10
    exploration:
      temperature: 300.0
    sampling:
      method: 'random'
      fraction: 0.8
    """
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    result = runner.invoke(app, ["--config", str(config_file)])

    assert result.exit_code == 0
    mock_pipeline_runner.assert_called_once()


def test_cli_file_not_found():
    """Test the CLI with a non-existent configuration file."""
    result = runner.invoke(app, ["--config", "non_existent_config.yaml"])

    assert result.exit_code == 1
    # Sanitize the output to handle rich formatting quirks
    sanitized_output = " ".join(result.stdout.strip().split())
    assert "Error" in sanitized_output
    assert "was not found" in sanitized_output
