"""End-to-end tests for the refactored Cycle 01 workflow."""

import pytest
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
import os
import yaml

from mlip_autopipec.application.cli import cli

SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =     -11.45567302 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000   -0.00000000

total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
      0.00000004   -0.00000000    0.00000000
      0.00000000    0.00000004    0.00000000
      0.00000000    0.00000000    0.00000004
"""

@pytest.fixture
def test_dir(tmp_path):
    config_data = {
        'db_path': str(tmp_path / 'test.db'),
        'model_path': str(tmp_path / 'model.pt')
    }
    with open(tmp_path / "config.yaml", "w") as f:
        yaml.dump(config_data, f)

    poscar_content = "Ar\n1.0\n10 0 0\n0 10 0\n0 0 10\n1\nDirect\n0 0 0"
    with open(tmp_path / "ar.poscar", "w") as f:
        f.write(poscar_content)

    return tmp_path

@patch('subprocess.run')
@patch('mlip_autopipec.application.services.torch.save')
def test_e2e_workflow_success(mock_torch_save, mock_subprocess_run, test_dir):
    mock_subprocess_run.return_value = MagicMock(stdout=SAMPLE_SUCCESS_OUTPUT, returncode=0)
    runner = CliRunner()
    config_path = str(test_dir / "config.yaml")

    result = runner.invoke(
        cli,
        ["run-cycle-01", "--config-file", config_path, "--input-dir", str(test_dir)],
    )

    assert result.exit_code == 0
    assert "--- MLIP-AutoPipe Cycle 01 Finished Successfully ---" in result.output
    mock_subprocess_run.assert_called_once()
    mock_torch_save.assert_called_once()
