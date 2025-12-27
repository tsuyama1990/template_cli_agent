"""Unit tests for the AseDB wrapper."""
import pytest
from ase.build import bulk
from pathlib import Path
import numpy as np

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


@pytest.fixture
def db_path(tmp_path: Path) -> Path:
    """Provides a temporary path for the test database."""
    return tmp_path / "test.db"


def test_asedb_write_and_get_successful(db_path: Path):
    """
    Tests writing and reading a successful DFT result.
    """
    # Arrange
    db = AseDB(db_path)
    si_atoms = bulk("Si", "diamond", a=5.43)
    # Correct 3x3 stress tensor
    stress_tensor = np.array([[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]])
    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        stress=stress_tensor.tolist(), # Store as list in Pydantic model
        was_successful=True,
    )

    # Act
    db_id = db.write(si_atoms, dft_result)
    retrieved_atoms, retrieved_result = db.get(db_id)

    # Assert
    assert db_id is not None
    assert retrieved_result.was_successful is True
    assert retrieved_result.error_message is None
    assert retrieved_result.total_energy_ev == -100.0

    # Check that forces and stress were retrieved correctly from the calculator
    assert np.allclose(retrieved_atoms.get_forces(), dft_result.forces)
    assert np.allclose(retrieved_atoms.get_stress(voigt=False), stress_tensor)

    # Verify the pydantic model also has the correct lists
    assert np.allclose(retrieved_result.forces, dft_result.forces)
    assert np.allclose(retrieved_result.stress, dft_result.stress)


def test_asedb_write_and_get_failed(db_path: Path):
    """
    Tests writing and reading a failed DFT result.
    """
    # Arrange
    db = AseDB(db_path)
    nacl_atoms = bulk("NaCl", "rocksalt", a=5.64)
    dft_result = DFTResult(
        was_successful=False,
        error_message="SCF failed to converge",
    )

    # Act
    db_id = db.write(nacl_atoms, dft_result)
    retrieved_atoms, retrieved_result = db.get(db_id)

    # Assert
    assert db_id is not None
    assert retrieved_result.was_successful is False
    assert retrieved_result.error_message == "SCF failed to converge"

    # For a failed run, the calculation results should be None
    assert retrieved_result.total_energy_ev is None
    assert retrieved_result.forces is None
    assert retrieved_result.stress is None

    # The Atoms object itself should not have a calculator attached
    assert not hasattr(retrieved_atoms, 'calc') or retrieved_atoms.calc is None
