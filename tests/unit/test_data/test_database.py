# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
from pathlib import Path

import numpy as np
import pytest
from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


@pytest.fixture
def temp_db(tmp_path: Path) -> Path:
    """Create a temporary database file for testing."""
    return tmp_path / "test.db"


def test_asedb_write_and_get_successful_result(temp_db: Path):
    """Verify writing and reading a successful DFT result."""
    db = AseDB(temp_db)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    dft_result = DFTResult(
        total_energy_ev=-31.0,
        forces=[[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]],
        stress=[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        was_successful=True,
    )

    db_id = db.write(atoms, dft_result)
    assert isinstance(db_id, int)

    retrieved_atoms, kvp = db.get(db_id)

    # Verify key-value pairs
    assert kvp["was_successful"] is True
    assert "error_message" not in kvp

    # Verify calculator data
    assert retrieved_atoms.calc is not None
    assert np.isclose(retrieved_atoms.get_potential_energy(), dft_result.total_energy_ev)
    assert np.allclose(retrieved_atoms.get_forces(), dft_result.forces)
    # ASE stores stress as a 6-element Voigt vector
    assert retrieved_atoms.get_stress().shape == (6,)


def test_asedb_write_and_get_failed_result(temp_db: Path):
    """Verify writing and reading a failed DFT result."""
    db = AseDB(temp_db)
    atoms = Atoms("O")
    dft_result = DFTResult(
        total_energy_ev=0.0,
        forces=[],
        stress=[],
        was_successful=False,
        error_message="SCF failed to converge",
    )

    db_id = db.write(atoms, dft_result)
    retrieved_atoms, kvp = db.get(db_id)

    # Verify key-value pairs for failure case
    assert kvp["was_successful"] is False
    assert kvp["error_message"] == "SCF failed to converge"

    # Verify no calculator is attached for a failed run
    assert retrieved_atoms.calc is None
