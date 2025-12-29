import os

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


@pytest.fixture
def temp_db(tmp_path):
    """Pytest fixture to create a temporary database for testing."""
    db_path = os.path.join(tmp_path, "test.db")
    yield AseDB(db_path)


def test_write_and_get_successful_result(temp_db):
    """
    Tests that a successful DFT result can be written to and retrieved from the database,
    with results attached to a SinglePointCalculator.
    """
    # 1. Setup
    atoms = bulk("Si", "diamond", a=5.43)
    stress_3x3 = [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
    successful_result = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0, 0.0, 0.0]] * 2,
        stress=stress_3x3,
        was_successful=True
    )

    # 2. Execute
    db_id = temp_db.write(atoms, successful_result)
    retrieved_atoms, retrieved_kvp = temp_db.get(db_id)

    # 3. Verify
    assert db_id == 1
    assert retrieved_kvp["was_successful"] is True
    assert "error_message" not in retrieved_kvp

    # Verify the calculator holds the results
    assert isinstance(retrieved_atoms.calc, SinglePointCalculator)
    assert retrieved_atoms.get_potential_energy() == successful_result.total_energy_ev
    assert np.allclose(retrieved_atoms.get_forces(), successful_result.forces)

    # ASE stores stress in Voigt form (6 elements), so convert the original for comparison
    expected_stress_voigt = full_3x3_to_voigt_6_stress(stress_3x3)
    assert np.allclose(retrieved_atoms.get_stress(), expected_stress_voigt)


def test_write_and_get_failed_result(temp_db):
    """
    Tests that a failed DFT result is correctly recorded in the database,
    with no calculator attached.
    """
    # 1. Setup
    atoms = bulk("Ar", "fcc", a=3.8)
    failed_result = DFTResult(
        total_energy_ev=0.0,
        forces=[],
        stress=[],
        was_successful=False,
        error_message="SCF failed to converge"
    )

    # 2. Execute
    db_id = temp_db.write(atoms, failed_result)
    retrieved_atoms, retrieved_kvp = temp_db.get(db_id)

    # 3. Verify
    assert db_id == 1
    assert retrieved_kvp["was_successful"] is False
    assert retrieved_kvp["error_message"] == "SCF failed to converge"

    # Ensure no calculator is attached on failure
    assert retrieved_atoms.calc is None
    # Check that energy/forces are not present in kvp either
    assert "total_energy_ev" not in retrieved_kvp


def test_get_nonexistent_id(temp_db):
    """
    Tests that querying a non-existent ID raises a KeyError with the correct message.
    """
    # 1. Setup - DB is empty

    # 2. Execute & Verify
    with pytest.raises(KeyError, match="No entry found with id 1"):
        temp_db.get(1)
