"""Unit tests for the CLI."""
from unittest.mock import patch

from typer.testing import CliRunner

from mlip_autopipec.cli import app

runner = CliRunner()


def test_cli_successful_invocation(tmp_path):
    """Test a successful invocation of the CLI."""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        """
project_name: uat_test_fept
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
    )

    with patch("mlip_autopipec.cli.PipelineRunner") as mock_runner:
        result = runner.invoke(app, ["--config", str(config_file)])
        assert result.exit_code == 0
        mock_runner.assert_called_once()
        mock_runner.return_value.run.assert_called_once()


def test_cli_file_not_found():
    """Test the CLI's behavior when the config file is not found."""
    result = runner.invoke(app, ["--config", "non_existent_config.yaml"])
    assert result.exit_code == 1

    # Sanitize the output to make the assertion more robust
    output_words = " ".join(result.output.strip().split())

    assert "Error" in output_words
    assert "not found" in output_words
