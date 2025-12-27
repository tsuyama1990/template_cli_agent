from pathlib import Path

from click.testing import CliRunner

from mlip_autopipec.main_cycle01 import main


def test_e2e_workflow():
    """Tests the full end-to-end workflow."""
    runner = CliRunner()
    config_path = Path("tests/unit/test_config.yaml")
    structure_path = Path("tests/e2e/test_structure.cif")
    db_path = Path("mlip.db")
    model_path = Path("models/test_model.pt")

    # Clean up previous runs
    if db_path.exists():
        db_path.unlink()
    if model_path.exists():
        model_path.unlink()

    result = runner.invoke(
        main,
        [
            "--config",
            str(config_path),
            "--structure",
            str(structure_path),
        ],
    )

    assert result.exit_code == 0
    assert "Workflow complete" in result.output
    assert db_path.exists()
    assert model_path.exists()

    # Clean up after the test
    db_path.unlink()
    model_path.unlink()
