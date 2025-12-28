# Description: End-to-end test for the Cycle 01 workflow.
import os
from pathlib import Path
from unittest.mock import patch

import pytest
from ase.io import write
from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import main


@pytest.fixture
def test_data(tmpdir):
    """Sets up the necessary files for an E2E test run."""
    # Create a simple structure file
    temp_path = Path(str(tmpdir))
    structure_path = temp_path / "Si.cif"
    from ase.build import bulk

    atoms = bulk("Si")
    write(str(structure_path), atoms)

    # Create a config file
    config_path = temp_path / "config.yaml"
    config_content = """
labelling_engine:
  qe_command: "mock_pw.sh -in QE_input.in"
training_engine:
  model_type: "None"
  learning_rate: 0.01
  num_epochs: 2
  r_cut: 4.0
  delta_learn: true
  baseline_potential: "ZBL"
"""
    config_path.write_text(config_content)

    # Create a mock QE executable and a mock output file
    mock_pw_path = temp_path / "mock_pw.sh"
    mock_output_content = """
     JOB DONE.
     !    total energy              =     -10.0 Ry
     Forces acting on atoms (Ry/au):
     atom    1 type  1   force =     0.0   0.0   0.0
     atom    2 type  1   force =     0.0   0.0   0.0
    """
    # The mock script just prints the mock output to stdout
    mock_pw_path.write_text(f"#!/bin/bash\necho \"{mock_output_content}\"")
    os.chmod(mock_pw_path, 0o755)

    # Add the temp directory to the PATH so the mock executable can be found
    original_path = os.environ["PATH"]
    os.environ["PATH"] = f"{str(temp_path)}:{original_path}"

    yield config_path, structure_path

    # Teardown: Restore the original PATH
    os.environ["PATH"] = original_path


def test_cycle01_e2e_workflow_success(test_data):
    """
    Tests the full end-to-end workflow for a successful run.
    """
    config_path, structure_path = test_data
    runner = CliRunner()

    # We need to patch the TrainingEngine's use of GPU if not available
    with patch(
        "mlip_autopipec.modules.d_training_engine.torch.cuda.is_available",
        return_value=False,
    ):
        result = runner.invoke(
            main,
            [
                "--config",
                str(config_path),
                "--structure",
                str(structure_path),
            ],
            catch_exceptions=False,
        )

    assert result.exit_code == 0
    assert "Workflow complete" in result.output
    assert Path("models/placeholder_model.txt").exists()
    assert Path("mlip.db").exists()

    # Cleanup created files
    os.remove("models/placeholder_model.txt")
    os.remove("mlip.db")
