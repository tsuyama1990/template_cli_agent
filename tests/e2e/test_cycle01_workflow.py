# tests/e2e/test_cycle01_workflow.py

import subprocess
from pathlib import Path

import pytest
import yaml
from ase import Atoms
from ase.io import write
from click.testing import CliRunner

from mlip_autopipec.cli import main

# Sample QE output for mocking the DFT calculation
SAMPLE_QE_OUTPUT = """
!    total energy              =     -15.85217439 Ry
Forces acting on atoms (cartesian axes, Ry/au):
     atom    1   force =     -0.00000014    -0.00000014    -0.00000014
total   stress  (Ry/bohr**3)                (kbar)     P=      -0.34
      -0.00000215    -0.00000000     0.00000000
      -0.00000000    -0.00000215     0.00000000
       0.00000000     0.00000000    -0.00000215
"""


@pytest.fixture
def test_workflow_files(tmp_path):
    """Fixture to create the necessary input files for an E2E test."""
    # 1. Create config.yaml
    config_data = {
        "system": {"elements": ["Si"]},
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "pw.x -v",
            "pseudopotentials": "SSSP",
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 1.0,
            "smearing": "mv",
            "degauss": 0.01,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": True,
            "base_potential": "lj_auto",
            "loss_weights": {"energy": 1.0, "force": 100.0},
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    # 2. Create initial_structures.xyz
    atoms = Atoms("Si", cell=[1, 1, 1], pbc=True)
    structures_path = tmp_path / "initial_structures.xyz"
    write(structures_path, atoms)

    # 3. Define database path
    db_path = tmp_path / "test.db"

    return config_path, structures_path, db_path


def test_cycle01_full_workflow_e2e(mocker, test_workflow_files):
    """
    An end-to-end test for the full Cycle 01 workflow.
    It uses the Click CLI runner and mocks the external DFT call.
    """
    config_path, structures_path, db_path = test_workflow_files

    # Mock the subprocess call to avoid running actual QE
    mock_subprocess_run = mocker.patch(
        "subprocess.run",
        return_value=subprocess.CompletedProcess(
            args=[], returncode=0, stdout=SAMPLE_QE_OUTPUT, stderr=""
        ),
    )

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--config-file",
            str(config_path),
            "--database-file",
            str(db_path),
            "--input-file",
            str(structures_path),
        ],
        catch_exceptions=False,  # This helps in debugging
    )

    assert result.exit_code == 0
    assert "Initializing Workflow" in result.output
    assert "Added 1 structures to the database" in result.output
    assert "Starting Labeling Engine" in result.output
    assert "Successfully labeled structure 1" in result.output
    assert "Starting Training Engine" in result.output
    assert "Saved trained model to model.ace" in result.output

    # Verify that the subprocess was called as expected
    mock_subprocess_run.assert_called_once()
    command_args = mock_subprocess_run.call_args[0][0]
    assert command_args[0] == "pw.x"
    assert command_args[1] == "-v"

    # Verify that the final model was created
    assert (Path.cwd() / "model.ace").exists()
    # Clean up the created model file
    (Path.cwd() / "model.ace").unlink()
