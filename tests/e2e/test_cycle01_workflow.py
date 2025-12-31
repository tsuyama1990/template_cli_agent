import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import run_cycle
from mlip_autopipec.data.database import AseDBWrapper


@pytest.fixture
def e2e_test_setup(tmp_path, mocker):
    # Mock the subprocess call to avoid running real DFT
    mocker.patch("subprocess.run")

    # Mock torch.save to avoid creating a real model file
    mocker.patch("torch.save")

    # Create a temporary config file
    db_path = tmp_path / "test_e2e.db"
    config = {
        "database_path": str(db_path),
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "pw.x",
            "pseudopotentials": "SSSP",
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 2.5,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": False,
            "loss_weights": {"energy": 1.0, "forces": 1.0},
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    # Create and populate a temporary database
    db = AseDBWrapper(str(db_path))
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], cell=[10, 10, 10])
    db.add_atoms(atoms, labeled=False)

    return config_path


def test_cycle01_e2e_workflow(e2e_test_setup):
    # Arrange
    config_path = e2e_test_setup
    runner = CliRunner()

    # Act
    result = runner.invoke(run_cycle, ["--config", str(config_path)])

    # Assert
    assert result.exit_code == 0
    assert "Cycle finished successfully" in result.output
    assert "Starting Labeling Engine" in result.output
    assert "Starting Training Engine" in result.output
