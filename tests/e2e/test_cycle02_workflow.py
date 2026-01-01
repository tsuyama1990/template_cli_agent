import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml
from click.testing import CliRunner
from ase import Atoms

from mlip_autopipec.cli import cli


@pytest.fixture
def runner() -> CliRunner:
    """Provides a Click test runner."""
    return CliRunner()


# To avoid running real DFT and training, we mock the execute methods
# of the engines and the database calls that provide data to them.
@patch("mlip_autopipec.data.database.AseDBWrapper.get_all_labeled_rows")
@patch("mlip_autopipec.data.database.AseDBWrapper.get_rows_to_label")
@patch("mlip_autopipec.modules.c_labeling_engine.LabelingEngine.execute")
@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine.execute")
def test_cycle02_full_workflow_e2e(
    mock_training_execute: MagicMock,
    mock_labeling_execute: MagicMock,
    mock_get_rows_to_label: MagicMock,
    mock_get_all_labeled_rows: MagicMock,
    runner: CliRunner,
    tmp_path: Path,
):
    """
    Tests the full 'Generate -> Label -> Train' workflow from the CLI.
    """
    # 1. Arrange:
    input_yaml_path = tmp_path / "input.yaml"
    input_data = {"system": {"elements": ["Cu", "Au"], "composition": "CuAu"}}
    with open(input_yaml_path, "w") as f:
        yaml.dump(input_data, f)

    # Configure the mock database methods to simulate a successful workflow
    # a) Provide one structure that needs labeling
    mock_row = MagicMock()
    mock_row.toatoms.return_value = Atoms("X")
    mock_get_rows_to_label.return_value = [mock_row]
    # b) Provide one labeled structure for the training engine
    mock_get_all_labeled_rows.return_value = [mock_row]

    os.chdir(tmp_path)

    # 2. Act:
    result = runner.invoke(
        cli, ["run", "--input", str(input_yaml_path)], catch_exceptions=False
    )

    # 3. Assert:
    assert result.exit_code == 0
    assert "Workflow finished successfully!" in result.output

    output_config_path = tmp_path / "exec_config_dump.yaml"
    db_path = tmp_path / "mlip_autopipec.db"
    assert output_config_path.exists()
    assert db_path.exists()

    with open(output_config_path) as f:
        expanded_config = yaml.safe_load(f)
    assert expanded_config["system"]["structure_type"] == "alloy"

    # Assert that the core modules were called.
    mock_labeling_execute.assert_called_once()
    mock_training_execute.assert_called_once()
