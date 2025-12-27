import torch.serialization
torch.serialization.add_safe_globals([slice])
import pytest
import os
import yaml
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from src.mlip_autopipec.main import main

# Sample successful Quantum Espresso output for the E2E test
SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =    -111.92138384 Ry
     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =    -0.00000002    -0.00000002    -0.00000002
     atom    2 type  1   force =     0.00000002     0.00000002     0.00000002

     Total Force =     0.00000004     Total SCF correction =     0.00000000

          total stress  (Ry/bohr**3)            (kbar)       P=   -2.29
   -0.00001452   0.00000000   0.00000000      -0.23      0.00      0.00
    0.00000000  -0.00001452   0.00000000       0.00     -0.23      0.00
    0.00000000   0.00000000  -0.00001452       0.00      0.00     -0.23
"""

@pytest.fixture
def test_files(tmp_path):
    """Creates dummy config and structure files for the E2E test."""
    config_dict = {
        "labelling_engine": {
            "qe_command": "pw.x",
            "pseudo_dir": "/dev/null",
            "pseudopotentials": {"Si": "Si.upf"},
        },
        "training_engine": {
            "model_type": "mace",
            "learning_rate": 0.01,
            "num_epochs": 1,
            "r_cut": 5.0,
            "delta_learn": True,
            "baseline_potential": "lennard_jones",
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)

    # Create a simple CIF file for a 2-atom Si cell
    cif_content = """
data_Si
_cell_length_a    5.43
_cell_length_b    5.43
_cell_length_c    5.43
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
"""
    structure_path = tmp_path / "si.cif"
    with open(structure_path, "w") as f:
        f.write(cif_content)

    return config_path, structure_path

@patch("subprocess.run")
def test_e2e_cycle01_workflow_success(mock_subprocess_run, test_files):
    """
    Tests the full end-to-end workflow for Cycle 01 using the CLI.
    Mocks the DFT calculation but uses the real training engine.
    """
    mock_subprocess_run.return_value = MagicMock(
        returncode=0, stdout=SAMPLE_SUCCESS_OUTPUT, stderr=""
    )

    config_path, structure_path = test_files

    runner = CliRunner()
    result = runner.invoke(
        main, ["--config", str(config_path), "--structure", str(structure_path)]
    )

    assert result.exit_code == 0
    assert "Workflow complete" in result.output

    # Check that a model file was actually created
    model_path = "trained_model.pt"
    assert os.path.exists(model_path)
    assert os.path.getsize(model_path) > 0

    # Cleanup the created files
    if os.path.exists("mlip.db"):
        os.remove("mlip.db")
    if os.path.exists(model_path):
        os.remove(model_path)
