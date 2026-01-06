# -*- coding: utf-8 -*-
"""End-to-end tests for the MLIP-AutoPipe pipeline."""
from __future__ import annotations

from pathlib import Path
import pytest
from typer.testing import CliRunner
from mlip_autopipec.cli import app

runner = CliRunner()


def test_full_pipeline_run(tmp_path: Path) -> None:
    """
    Test a full, end-to-end run of the minimal pipeline.

    This test runs the CLI with a valid configuration and verifies that the
    final database file is created, confirming that all components integrate
    correctly.
    """
    # 1. Create a valid configuration file in a temporary directory
    config_content = """
    project_name: integration_test
    system:
      elements: ['Ni', 'Al']
      composition: {'Ni': 0.8, 'Al': 0.2}
      lattice: 'fcc'
      num_structures: 2
    exploration:
      temperature: 500.0
    sampling:
      method: 'random'
      fraction: 0.5
    """

    with runner.isolated_filesystem(temp_dir=tmp_path) as td:
        config_file_iso = Path(td) / "config.yaml"
        config_file_iso.write_text(config_content)

        # 2. Run the CLI command, catching the expected SystemExit
        with pytest.raises(SystemExit) as e:
            runner.invoke(app, ["run", "--config", str(config_file_iso)], catch_exceptions=False)

        # 3. Assert that the exit code is 0 (success)
        assert e.value.code == 0, f"CLI exited with a non-zero exit code: {e.value.code}"

        # 4. Verify that the output database was created
        output_db_iso = Path(td) / "integration_test.db"
        assert output_db_iso.exists(), "Output database was not created."
