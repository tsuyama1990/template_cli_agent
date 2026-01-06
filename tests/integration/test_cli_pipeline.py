# tests/integration/test_cli_pipeline.py
import tempfile
from pathlib import Path
from typing import Generator

import pytest
import yaml
from ase.db import connect

from mlip_autopipec.cli.main import run_pipeline


@pytest.fixture
def temp_config_and_db() -> Generator[tuple[Path, Path], None, None]:
    """
    Creates a temporary directory, a valid YAML config file, and a path for the
    output database. Cleans up afterward.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = Path(tmpdir) / "config.yaml"
        db_path = Path(tmpdir) / "test_output.db"

        config_data = {
            "system": {
                "elements": ["Cu", "Au"],
                "composition": {"Cu": 0.5, "Au": 0.5},
                "supercell_size": [1, 1, 1],
            },
            "exploration": {
                "temperature_k": 10.0,  # Low temp for stability
                "pressure_gpa": 0,
                "timestep_fs": 0.5,
                "n_steps": 5,  # Very short run
            },
            "sampling": {"method": "random", "n_samples": 3},
            "db_path": str(db_path),
        }

        with config_path.open("w") as f:
            yaml.dump(config_data, f)

        yield config_path, db_path

        # Cleanup is handled by TemporaryDirectory context manager


def test_cli_pipeline_end_to_end(temp_config_and_db: tuple[Path, Path]) -> None:
    """
    Tests the full CLI pipeline from a config file to a final database.
    """
    config_path, db_path = temp_config_and_db

    # Invoke the command function directly
    run_pipeline(str(config_path))

    # Verify the output
    assert db_path.exists(), "Database file was not created"

    # Check the database content
    db = connect(db_path)  # type: ignore[no-untyped-call]
    # The number of structures in the DB should match n_samples
    assert len(db) == 3
    # Check that one of the rows has the correct elements
    row = db.get(id=1)
    assert "Cu" in row.symbols or "Au" in row.symbols
