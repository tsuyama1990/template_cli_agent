# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
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


@patch("mlip_autopipec.modules.d_training_engine.torch.save")
@patch("mlip_autopipec.modules.d_training_engine.MACE")
@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
def test_e2e_successful_workflow(
    mock_subprocess_run: MagicMock,
    mock_mace: MagicMock,
    mock_torch_save: MagicMock,
    test_files: dict,
):
    """
    Tests the full end-to-end workflow for a successful run (UAT-C01-01, UAT-C01-02),
    mocking the external QE call and the MACE model itself to bypass library bugs.
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

    # Mock the MACE model to avoid the library's internal TypeError
    mock_model_instance = MagicMock()
    mock_model_instance.parameters.return_value = [
        torch.nn.Parameter(torch.randn(1, requires_grad=True))
    ]
    mock_model_instance.return_value = {
        "energy": torch.tensor([0.0], requires_grad=True, dtype=torch.float64),
        "forces": torch.tensor([[0.0, 0.0, 0.0]], requires_grad=True, dtype=torch.float64),
    }
    mock_mace.return_value = mock_model_instance

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

    # UAT-C01-01: Check for model file
    mock_torch_save.assert_called_once()
    model_path_match = re.search(r"Model saved to: (models/trained_model.pt)", result.output)
    assert model_path_match is not None, "Model path not found in output."

    # UAT-C01-02: Check database contents
    db_file = test_files["db"]
    assert db_file.exists()
    with connect(db_file) as db:
        assert len(db) == 1
        row = db.get(id=1)
        atoms = row.toatoms()
        assert np.allclose(atoms.get_forces(), 0.0, atol=1e-5)
        assert row.energy < 0


@patch("mlip_autopipec.modules.c_labelling_engine.subprocess.run")
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

    # UAT-C01-04: Check for no model and correct DB entry
    assert "Model saved to:" not in result.output
    assert "ERROR: Labelling failed" in result.output

    db_file = test_files["db"]
    assert db_file.exists()
    with connect(db_file) as db:
        assert len(db) == 1
        row = db.get(id=1)
        assert row.key_value_pairs.get("was_successful") is False
        assert "SCF failed to converge" in row.key_value_pairs.get("error_message", "")
