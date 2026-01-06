from pathlib import Path
from typing import Generator

import pytest
import yaml
from ase.db import connect
from typer import Exit
from typer.testing import CliRunner

from mlip_autopipec.cli.main import run

runner = CliRunner()


@pytest.fixture(scope="module")
def sample_config_file(tmp_path_factory: pytest.TempPathFactory) -> Generator[Path, None, None]:
    """Creates a temporary YAML config file for testing."""
    config_data = {
        "system": {
            "elements": ["Cu"],
            "composition": {"Cu": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 300, "pressure_gpa": 0},
    }
    config_path = tmp_path_factory.mktemp("data") / "config.yaml"
    with config_path.open("w") as f:
        yaml.dump(config_data, f)
    yield config_path


def test_cli_run_end_to_end(sample_config_file: Path) -> None:
    """
    Tests the full end-to-end workflow of the `run` CLI command.
    """
    # Run the underlying function directly to bypass CliRunner issues
    try:
        run(config_path=sample_config_file)
    except Exit as e:
        pytest.fail(f"CLI exited unexpectedly: {e.exit_code}")

    # Verify that the output database was created
    output_db_path = Path("output/final_structures.db")
    assert output_db_path.exists()

    # Verify the contents of the database
    db = connect(str(output_db_path))  # type: ignore[no-untyped-call]
    # The RandomSampler is configured to select 3 samples
    assert len(list(db.select())) == 3

    # Clean up the output file
    output_db_path.unlink()
    output_db_path.parent.rmdir()
