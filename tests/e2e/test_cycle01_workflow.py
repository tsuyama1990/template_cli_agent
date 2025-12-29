# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
import hashlib
import re
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import torch
from ase.db import connect
from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import main

# Mark all tests in this file as E2E
pytestmark = pytest.mark.e2e


@pytest.fixture
def test_files(tmp_path: Path) -> dict:
    """Sets up common files needed for E2E tests."""
    config_content = """
qe_command: "pw.x"
pseudo_potentials:
  Si: "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
dft_params:
  ecutwfc: 30.0
  k_points: [1, 1, 1]
training:
  model_type: "mace"
  learning_rate: 0.01
  num_epochs: 2
  r_cut: 4.0
  delta_learn: false
  baseline_potential: "lj"
"""
    cif_content = """
data_Si
_cell_length_a   5.43
_cell_length_b   5.43
_cell_length_c   5.43
_cell_angle_alpha   90
_cell_angle_beta   90
_cell_angle_gamma   90
loop_
  _atom_site_label
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  Si   0.00000   0.00000   0.00000
  Si   0.25000   0.25000   0.25000
"""
    pseudo_content = "dummy pseudo content"

    config_file = tmp_path / "test_config.yaml"
    cif_file = tmp_path / "si_bulk.cif"
    pseudo_file = tmp_path / "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
    db_file = tmp_path / "test.db"

    config_file.write_text(config_content)
    cif_file.write_text(cif_content)
    pseudo_file.write_text(pseudo_content)

    return {"config": config_file, "cif": cif_file, "db": db_file}


def file_md5(path: Path) -> str:
    """Computes the MD5 hash of a file."""
    if not path.exists():
        return ""
    return hashlib.md5(path.read_bytes()).hexdigest()


class DummyMACE(torch.nn.Module):
    """A minimal, serializable torch module to stand in for the real MACE model."""

    def __init__(self, **kwargs):
        super().__init__()
        self.linear = torch.nn.Linear(1, 1)

    def forward(self, data: dict) -> dict:
        # A dummy calculation that uses the input to produce some output
        # Ensure output tensors require gradients for the loss function
        num_atoms = len(data.get("ptr", [0, 1])) - 1
        return {
            "energy": torch.tensor([0.5] * num_atoms, requires_grad=True, dtype=torch.float64),
            "forces": torch.zeros(
                (num_atoms, 3), requires_grad=True, dtype=torch.float64
            ),
        }


@patch("mlip_autopipec.infrastructure.process_runner.subprocess.run")
def test_e2e_successful_workflow(mock_subprocess_run: MagicMock, test_files: dict):
    """
    Tests the full end-to-end workflow for a successful run (UAT-C01-01, UAT-C01-02).
    """
    mock_subprocess_run.return_value.stdout = """
!    total energy              =     -78.0 Ry
Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =     0.0000001  0.0  0.0
     atom    2 type  1   force =    -0.0000001  0.0  0.0
total stress  (Ry/bohr**3)                (kbar)     P=       0.0
      0.0   0.0   0.0    0.0   0.0   0.0
    """
    mock_subprocess_run.return_value.returncode = 0

    with patch(
        "mlip_autopipec.modules.d_training_engine.MACE", new=DummyMACE
    ) as mock_mace:
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "--config",
                str(test_files["config"]),
                "--structure",
                str(test_files["cif"]),
                "--db",
                str(test_files["db"]),
            ],
            catch_exceptions=False,
        )

    assert result.exit_code == 0, f"CLI crashed with output:\\n{result.output}"

    model_path_match = re.search(
        r"Model saved to: (models/trained_model.pt)", result.output
    )
    assert model_path_match is not None, "Model path not found in output."

    db_file = test_files["db"]
    assert db_file.exists()
    with connect(db_file) as db:
        assert len(db) == 1
        row = db.get(id=1)
        atoms = row.toatoms()
        assert np.allclose(atoms.get_forces(), 0.0, atol=1e-5)
        assert row.energy < 0


@patch("mlip_autopipec.infrastructure.process_runner.subprocess.run")
def test_e2e_dft_failure_workflow(mock_subprocess_run: MagicMock, test_files: dict):
    """
    Tests the workflow's handling of a failed DFT calculation (UAT-C01-04).
    """
    mock_subprocess_run.return_value.stdout = "SCF NOT CONVERGED"
    mock_subprocess_run.return_value.returncode = 0

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--config",
            str(test_files["config"]),
            "--structure",
            str(test_files["cif"]),
            "--db",
            str(test_files["db"]),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, f"CLI crashed with output:\\n{result.output}"
    assert "Model saved to:" not in result.output
    assert "ERROR: Labelling failed" in result.output

    db_file = test_files["db"]
    assert db_file.exists()
    with connect(db_file) as db:
        assert len(db) == 1
        row = db.get(id=1)
        assert row.key_value_pairs.get("was_successful") is False
        assert "SCF failed to converge" in row.key_value_pairs.get("error_message", "")


@patch("mlip_autopipec.infrastructure.process_runner.subprocess.run")
def test_e2e_delta_learning_toggle(
    mock_subprocess_run: MagicMock, test_files: dict, tmp_path: Path
):
    """
    Tests that the delta_learn flag correctly alters the training process (UAT-C01-03).
    """
    mock_subprocess_run.return_value.stdout = """
!    total energy              =     -78.0 Ry
Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =     0.1  0.1  0.1
     atom    2 type  1   force =    -0.1 -0.1 -0.1
total stress  (Ry/bohr**3)                (kbar)     P=       0.0
      0.0   0.0   0.0    0.0   0.0   0.0
    """
    mock_subprocess_run.return_value.returncode = 0

    with patch(
        "mlip_autopipec.modules.d_training_engine.MACE", new=DummyMACE
    ) as mock_mace:
        # Run 1: Delta Learning Disabled
        runner = CliRunner()
        result_direct = runner.invoke(
            main,
            [
                "--config",
                str(test_files["config"]),
                "--structure",
                str(test_files["cif"]),
                "--db",
                str(test_files["db"]),
            ],
            catch_exceptions=False,
        )
        assert result_direct.exit_code == 0
        model_path_direct = Path("models/trained_model.pt")
        assert model_path_direct.exists()
        hash_direct = file_md5(model_path_direct)

        # Run 2: Delta Learning Enabled
        config_delta_content = test_files["config"].read_text().replace("delta_learn: false", "delta_learn: true")
        config_file_delta = tmp_path / "test_config_delta.yaml"
        config_file_delta.write_text(config_delta_content)

        runner_delta = CliRunner()
        result_delta = runner_delta.invoke(
            main,
            [
                "--config",
                str(config_file_delta),
                "--structure",
                str(test_files["cif"]),
                "--db",
                str(test_files["db"]),
            ],
            catch_exceptions=False,
        )
        assert result_delta.exit_code == 0
        model_path_delta = Path("models/trained_model.pt")
        assert model_path_delta.exists()
        hash_delta = file_md5(model_path_delta)

    assert hash_direct != hash_delta, "Delta learning flag did not change the model output."
