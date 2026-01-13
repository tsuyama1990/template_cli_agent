from unittest.mock import patch

from typer.testing import CliRunner

from mlip_autopipec.cli import app

runner = CliRunner()


@patch("mlip_autopipec.cli.PipelineRunner")
def test_cli_successful_invocation(mock_pipeline_runner):
    """Tests that the CLI can be invoked successfully."""
    result = runner.invoke(app, ["--config", "tests/test_config.yaml"])
    assert result.exit_code == 0
    mock_pipeline_runner.assert_called_once()
    mock_pipeline_runner.return_value.run.assert_called_once()


def test_cli_file_not_found():
    """Tests that the CLI exits gracefully if the config file is not found."""
    result = runner.invoke(app, ["--config", "non_existent_file.yaml"])
    assert result.exit_code != 0
