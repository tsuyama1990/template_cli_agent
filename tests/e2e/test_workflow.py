"""End-to-end test for the Cycle 01 workflow."""
import pytest
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
from pathlib import Path
import subprocess
import yaml
from ase import Atoms

from mlip_autopipec.main import main

# Re-use the sample successful QE output for the mock
SAMPLE_QE_OUTPUT = """
!    total energy              =      -2.27598151 Ry
Forces acting on atoms (Ry/au):
     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000    0.00000000
total stress  (Ry/bohr**3)                (kbar)     P=        0.00
   0.00000000    0.00000000    0.00000000    0.00000000    0.00000000    0.00000000
JOB DONE.
"""

@pytest.fixture
def test_environment(tmp_path: Path) -> tuple[Path, Path]:
    """Creates a temporary environment with a config file and a structure file."""

    # Create a minimal config file
    config_data = {
        "qe_command": "pw.x",
        "pseudo_dir": "pseudos/",
        "pseudos": {"H": "H.UPF"},
        "kpts": [1, 1, 1],
        "ecutwfc": 30.0,
        "training": {
            "model_type": "mace",
            "learning_rate": 0.01,
            "num_epochs": 2, # Minimal epochs for a fast test
            "r_cut": 4.0,
            "delta_learn": True,
            "baseline_potential": "lj",
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    # Create a simple structure file
    h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure_path = tmp_path / "h2.xyz"
    h2.write(structure_path, format="extxyz")

    return config_path, structure_path


@patch("subprocess.run")
def test_e2e_workflow_success(mock_subprocess_run: MagicMock, test_environment: tuple[Path, Path]):
    """
    Tests the full end-to-end workflow from the CLI.
    Mocks the external Quantum Espresso call.
    """
    # Arrange
    config_path, structure_path = test_environment

    # Mock the subprocess call to return a successful QE run
    mock_subprocess_run.return_value = subprocess.CompletedProcess(
        args=["pw.x", "-in", "qe.in"],
        returncode=0,
        stdout=SAMPLE_QE_OUTPUT,
        stderr="",
    )

    runner = CliRunner()

    # Act
    with runner.isolated_filesystem() as fs:
        # We need to copy the files into the isolated filesystem
        isolated_config_path = Path(fs) / "config.yaml"
        isolated_structure_path = Path(fs) / "h2.xyz"
        config_path.rename(isolated_config_path)
        structure_path.rename(isolated_structure_path)

        result = runner.invoke(
            main,
            ["--config", str(isolated_config_path), "--structure", str(isolated_structure_path)],
        )

        # Assert
        assert result.exit_code == 0
        assert "Workflow Complete" in result.output
        assert "Trained model saved to" in result.output

        # Check that the model file was actually created
        model_file = Path(fs) / "trained_model.pt"
        assert model_file.exists()
        assert model_file.stat().st_size > 0

        # Check that the database was created
        db_file = Path(fs) / "mlip_pipe.db"
        assert db_file.exists()
