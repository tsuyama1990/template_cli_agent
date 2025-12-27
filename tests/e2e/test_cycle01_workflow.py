from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml
from ase.build import bulk
from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import main as cli_main

# Re-using the sample output from the unit test
SAMPLE_QE_SUCCESS_OUTPUT = """
     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =    -0.00000000    -0.00000000    -0.00000000
     atom    2 type  1   force =     0.00000000     0.00000000     0.00000000
!    total energy              =     -11.43845671 Ry
     total stress  (Ry/bohr**3)                (kbar)     P=  -13.71
     -0.00086555   0.00000000   0.00000000        -12.72    0.00    0.00
      0.00000000  -0.00086555   0.00000000          0.00  -12.72    0.00
      0.00000000   0.00000000  -0.00086555          0.00    0.00  -12.72
     convergence has been achieved in  8 iterations
     JOB DONE.
"""


@pytest.fixture
def test_setup(tmp_path: Path):
    """Creates a temporary directory with a config file, structure file, and mock DB."""
    db_path = tmp_path / "test.db"
    config_data = {
        "database": {"path": str(db_path)},
        "dft": {
            "qe_command": "pw.x",
            "parameters": {"system": {"ecutwfc": 60}},
            "pseudopotentials": {"Si": "Si.upf"},
        },
        "training": {
            "model_type": "mace",
            "learning_rate": 0.001,
            "num_epochs": 1,
            "r_cut": 4.5,
            "delta_learn": True,
            "baseline_potential": "lj",
        },
    }
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    structure_file = tmp_path / "si.cif"
    atoms = bulk("Si", "diamond", a=5.43)
    from ase.io import write

    write(str(structure_file), atoms)

    mock_db = MagicMock()
    mock_db.get.return_value = (
        atoms,
        {
            "total_energy_ev": -150.0,
            "forces": [[0, 0, 0]] * 2,
            "stress": [[0] * 3] * 3,
            "was_successful": True,
        },
    )

    return config_file, structure_file, mock_db, db_path


@patch(
    "mlip_autopipec.modules.d_training_engine.TrainingEngine.execute",
    return_value="models/trained_model.pt",
)
@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine._get_z_table")
@patch("mlip_autopipec.orchestrator_cycle01.AseDB")
@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
def test_e2e_cycle01_workflow_success(
    mock_subprocess_run, mock_asedb, mock_get_z_table, mock_training_execute, test_setup
):
    """
    Tests the full E2E workflow for a successful run using the CLI.
    Mocks the external subprocess call and the database connection.
    """
    config_file, structure_file, mock_db_instance, db_path = test_setup
    mock_asedb.return_value = mock_db_instance
    from mace.tools import AtomicNumberTable

    mock_get_z_table.return_value = AtomicNumberTable([14])  # Mock for Silicon

    # Arrange: Mock the side effects
    def mock_run(*args, **kwargs):
        cmd = args[0]
        output_path_str = cmd.split(" > ")[1]
        Path(output_path_str).write_text(SAMPLE_QE_SUCCESS_OUTPUT)
        return MagicMock(returncode=0)

    mock_subprocess_run.side_effect = mock_run

    # Act: Run the CLI
    runner = CliRunner()
    result = runner.invoke(
        cli_main,
        ["--config", str(config_file), "--structure", str(structure_file)],
    )

    # Assert
    assert result.exit_code == 0, result.output
    assert "Workflow complete." in result.output
