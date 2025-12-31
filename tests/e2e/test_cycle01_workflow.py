"""End-to-end tests for the Cycle 01 workflow."""

from unittest.mock import MagicMock, patch

import pytest
import yaml
from click.testing import CliRunner

from mlip_autopipec.main import cli

# Using a sample output that is known to be valid for the parser
SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =     -11.45567302 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000    0.00000000
     atom    2   force =     0.00000101    0.00000000    0.00000000

total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
      0.00000004   -0.00000000    0.00000000
      0.00000000    0.00000004    0.00000000
      0.00000000    0.00000000    0.00000004
"""


@pytest.fixture
def cycle01_test_dir(tmp_path):
    """Fixture to create a temporary directory with files for a Cycle 01 run."""
    config_data = {
        "db_path": str(tmp_path / "test.db"),
        "dft_command": "pw.x",
        "model_path": str(tmp_path / "model.pt"),
    }
    with open(tmp_path / "config.yaml", "w") as f:
        yaml.dump(config_data, f)

    poscar_content = """Ar2
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
Ar
2
Direct
0.0 0.0 0.0
0.0 0.0 0.1
"""
    with open(tmp_path / "ar2.poscar", "w") as f:
        f.write(poscar_content)

    return tmp_path


@patch("subprocess.run")
@patch("mlip_autopipec.modules.training_engine.MACE")
@patch("mlip_autopipec.modules.training_engine.torch.save")
def test_cycle01_e2e_workflow_success(
    mock_torch_save, mock_mace, mock_subprocess_run, cycle01_test_dir
):
    """
    Tests the full end-to-end workflow for Cycle 01 with a mocked DFT call.
    """
    # Arrange
    mock_subprocess_run.return_value = MagicMock(stdout=SAMPLE_SUCCESS_OUTPUT, returncode=0)
    mock_mace.return_value.state_dict.return_value = {"param": 1}

    runner = CliRunner()
    config_path = str(cycle01_test_dir / "config.yaml")
    model_path = str(cycle01_test_dir / "model.pt")

    # Act
    result = runner.invoke(
        cli,
        ["run-cycle-01", "--config-file", config_path, "--input-dir", str(cycle01_test_dir)],
        catch_exceptions=False,  # Crucial for debugging
    )

    # Assert
    assert result.exit_code == 0
    assert "--- MLIP-AutoPipe Cycle 01 Finished Successfully ---" in result.output
    assert "Loaded 1 structures into the database." in result.output
    assert "--- Labelling Engine Finished ---" in result.output
    assert "--- Training Engine Finished ---" in result.output

    mock_subprocess_run.assert_called_once()
    mock_torch_save.assert_called_once()
    assert mock_torch_save.call_args[0][1] == model_path
