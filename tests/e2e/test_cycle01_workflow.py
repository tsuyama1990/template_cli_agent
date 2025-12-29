from unittest.mock import Mock, patch

import numpy as np
import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.main_cycle01 import main as cli_main

# --- Test Data ---

# A simple, valid QE output for an H2 molecule
MOCK_QE_OUTPUT = """
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000   -0.00000000

     !    total energy              =     -2.30158634 Ry

     total stress  (Ry/bohr**3)  (kbar)     P=      -0.01
       -0.00010   0.00000   0.00000
        0.00000  -0.00010   0.00000
        0.00000   0.00000  -0.00010

     JOB DONE.
"""

# --- Fixtures ---

@pytest.fixture(scope="module")
def runner():
    """Fixture for invoking command-line interfaces."""
    return CliRunner()

@pytest.fixture
def test_files(tmp_path):
    """Fixture to create necessary files for an E2E test run in a temporary directory."""
    # Create structure file
    h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure_file = tmp_path / "h2.xyz"
    h2.write(structure_file)

    # Create config file
    config_data = {
        "qe_command": "pw.x", # This will be mocked
        "db_path": str(tmp_path / "test.db"),
        "dft_parameters": {
            "pseudopotentials": {"H": "H.UPF"},
            "k_points": [1, 1, 1],
            "ecutwfc": 30,
            "ecutrho": 120,
        },
        "training": {
            "model_type": "MACE",
            "learning_rate": 0.01,
            "num_epochs": 2,
            "r_cut": 3.0,
            "delta_learn": True,
            "baseline_potential": "LJ",
        },
    }
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    return config_file, structure_file, tmp_path

# --- E2E Test ---

@patch("subprocess.run")
def test_cycle01_e2e_successful_workflow_with_placeholder(mock_subprocess_run, runner, test_files):
    """
    Tests the full end-to-end workflow for Cycle 01 using the placeholder TrainingEngine.
    """
    # 1. Setup
    config_file, structure_file, temp_dir = test_files

    mock_process = Mock()
    mock_process.stdout = MOCK_QE_OUTPUT
    mock_process.stderr = ""
    mock_process.returncode = 0
    mock_subprocess_run.return_value = mock_process

    # 2. Execute
    result = runner.invoke(
        cli_main, ["--config", str(config_file), "--structure", str(structure_file)]
    )

    # 3. Verify
    assert result.exit_code == 0
    assert "Workflow complete" in result.output

    # Check for the placeholder model file
    model_path = temp_dir / "models" / "mace_model_placeholder.pt"
    assert model_path.exists()
    assert model_path.stat().st_size > 0

    db_path = temp_dir / "test.db"
    assert db_path.exists()
    db = AseDB(str(db_path))
    retrieved_atoms, kvp = db.get(1)

    assert kvp["was_successful"] is True
    assert pytest.approx(retrieved_atoms.get_potential_energy(), 0.1) == -31.3
    assert np.allclose(retrieved_atoms.get_forces(), 0.0, atol=1e-5)
