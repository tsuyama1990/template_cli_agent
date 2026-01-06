from pathlib import Path
from unittest.mock import patch

import pytest
import typer
from typer.testing import CliRunner

from mlip_autopipec.cli import app

runner = CliRunner()


def test_cli_successful_run(tmp_path: Path) -> None:
    """Test a successful run of the CLI."""
    config_content = """
project_name: test_cli
system:
  elements: [Fe, Pt]
  composition: {Fe: 0.5, Pt: 0.5}
  lattice: fcc
  num_structures: 1
exploration:
  temperature: 300
sampling:
  method: random
  fraction: 0.5
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    with patch("mlip_autopipec.cli.PipelineRunner") as mock_runner:
        result = runner.invoke(app, ["--config", str(config_file)])
        assert result.exit_code == 0
        mock_runner.assert_called_once()
        mock_runner.return_value.run.assert_called_once()
        assert "Pipeline finished successfully" in result.output


def test_cli_file_not_found() -> None:
    """Test the CLI with a non-existent config file."""
    result = runner.invoke(app, ["--config", "non_existent_file.yaml"])
    assert result.exit_code == 1
    assert "Error: Config file not found" in result.output


def test_cli_general_exception(tmp_path: Path) -> None:
    """Test that the CLI handles general exceptions gracefully."""
    config_content = """
project_name: test_cli
system:
  elements: [Fe, Pt]
  composition: {Fe: 0.5, Pt: 0.5}
  lattice: fcc
  num_structures: 1
exploration:
  temperature: 300
sampling:
  method: random
  fraction: 0.5
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    with patch("mlip_autopipec.cli.PipelineRunner") as mock_runner:
        mock_runner.return_value.run.side_effect = ValueError("Test exception")
        result = runner.invoke(app, ["--config", str(config_file)])
        assert result.exit_code == 1
        assert "An error occurred: Test exception" in result.output
