import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import cli
from mlip_autopipec.data.models import FullConfig


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def full_config_dict() -> dict:
    """Provides a valid, fully expanded config as a dictionary."""
    return {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1},
            "structure_type": "covalent",
        },
        "generation": {
            "generation_strategy": "random",
            "supercell_size": 8,
            "strains": [0.0],
        },
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "pw.x",
            "pseudopotentials": {"Si": "Si.UPF"},
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 3.0,
        },
        "mlip_training": {
            "model_type": "mace",
            "r_cut": 5.0,
            "delta_learning": False,
            "loss_weights": {"energy": 1.0, "forces": 1.0},
        },
    }


@patch("mlip_autopipec.data.database.AseDBWrapper.get_all_labeled_rows")
@patch("mlip_autopipec.data.database.AseDBWrapper.get_rows_to_label")
@patch("mlip_autopipec.services.config_expander.ConfigExpander.expand_config")
@patch("mlip_autopipec.modules.a_structure_generator.StructureGenerator.execute")
@patch("mlip_autopipec.modules.c_labeling_engine.LabelingEngine.execute")
@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine.execute")
def test_cycle01_logic_within_new_cli(
    mock_training_execute: MagicMock,
    mock_labeling_execute: MagicMock,
    mock_structure_generator_execute: MagicMock,
    mock_expand_config: MagicMock,
    mock_get_rows_to_label: MagicMock,
    mock_get_all_labeled_rows: MagicMock,
    runner: CliRunner,
    tmp_path: Path,
    full_config_dict: dict,
):
    """
    Tests that the core 'Label -> Train' logic from Cycle 1 still works
    within the new Cycle 2 CLI and orchestrator.
    """
    # Arrange
    mock_expand_config.return_value = FullConfig(**full_config_dict)

    # Simulate a successful workflow by mocking the database calls
    mock_row = MagicMock()
    mock_row.toatoms.return_value = Atoms("X")
    mock_get_rows_to_label.return_value = [mock_row]
    mock_get_all_labeled_rows.return_value = [mock_row]

    input_yaml_path = tmp_path / "input.yaml"
    input_yaml_path.touch()
    os.chdir(tmp_path)

    # Act
    result = runner.invoke(
        cli, ["run", "--input", str(input_yaml_path)], catch_exceptions=False
    )

    # Assert
    assert result.exit_code == 0
    assert "Workflow finished successfully!" in result.output
    mock_structure_generator_execute.assert_called_once()
    mock_labeling_execute.assert_called_once()
    mock_training_execute.assert_called_once()


def test_cli_run_file_not_found(runner: CliRunner):
    """Tests the CLI with a non-existent input file."""
    result = runner.invoke(
        cli, ["run", "--input", "nonexistent.yaml"], catch_exceptions=False
    )
    assert result.exit_code != 0
    assert "File 'nonexistent.yaml' does not exist." in result.output
