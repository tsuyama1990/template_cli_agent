from pathlib import Path

import pytest
from typer.testing import CliRunner

# from mlip_autopipec.cli.main import app

runner = CliRunner()

@pytest.fixture
def test_config_yaml(tmp_path: Path) -> Path:
    """Creates a temporary YAML config file for testing."""
    config_content = """
system:
  elements: [Si]
  composition: {Si: 1.0}
  supercell_size: [2, 2, 2]

exploration:
  temperature_k: 300.0
  pressure_gpa: 0.0
  time_step_fs: 1.0
  total_steps: 10 # Keep low for testing

sampling:
  method: random
  num_samples: 2
"""
    config_file = tmp_path / "test_config.yaml"
    config_file.write_text(config_content)
    return config_file

@pytest.mark.skip(reason="CLI implementation not yet available")
def test_cli_end_to_end_run(test_config_yaml: Path, tmp_path: Path):
    """
    Tests the full end-to-end pipeline via a CLI command.
    """
    # output_db = tmp_path / "results.db"
    #
    # result = runner.invoke(
    #     app,
    #     ["run", str(test_config_yaml), f"--db-path={output_db}"],
    #     catch_exceptions=False
    # )
    #
    # assert result.exit_code == 0
    # assert "Pipeline completed successfully" in result.output
    #
    # # Verify that the final database was created
    # assert output_db.exists()
    #
    # # Further assertions would connect to the DB and check its contents
    # # from ase.db import connect
    # # db = connect(output_db)
    # # assert len(db) == 2 # Check if the correct number of samples were saved
