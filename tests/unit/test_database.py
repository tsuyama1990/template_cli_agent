import pytest
import os
import numpy as np
from ase.atoms import Atoms
from ase.build import bulk
from src.mlip_autopipec.data.database import AseDB
from src.mlip_autopipec.data.models import DFTResult

@pytest.fixture
def temp_db():
    """Provides a temporary database path and ensures cleanup."""
    db_path = "test_temp.db"
    yield db_path
    if os.path.exists(db_path):
        os.remove(db_path)

def test_write_and_read_successful_result(temp_db):
    """
    Tests writing a successful DFT result to the database and reading it back.
    """
    db = AseDB(db_path=temp_db)
    atoms = bulk("Si", "diamond", a=5.43)

    dft_result = DFTResult(
        total_energy_ev=-100.0,
        forces=[[0.0, 0.0, 0.0]] * len(atoms),
        stress=[[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
        was_successful=True
    )

    # Write to DB
    db_id = db.write(atoms, dft_result)
    assert isinstance(db_id, int)

    # Read from DB
    retrieved_atoms, kvp = db.get(db_id)

    # Verify Atoms object
    assert np.array_equal(atoms.get_positions(), retrieved_atoms.get_positions())
    assert np.array_equal(atoms.get_cell(), retrieved_atoms.get_cell())
    assert atoms.get_chemical_symbols() == retrieved_atoms.get_chemical_symbols()

    # Verify calculator results attached to Atoms
    assert retrieved_atoms.get_potential_energy() == pytest.approx(dft_result.total_energy_ev)
    assert np.allclose(retrieved_atoms.get_forces(), dft_result.forces)
    # ASE stores stress as a 6-element Voigt vector, so we compare with the full 3x3 matrix
    assert np.allclose(retrieved_atoms.get_stress(voigt=False), dft_result.stress)

    # Verify key-value pairs
    assert kvp["total_energy_ev"] == pytest.approx(dft_result.total_energy_ev)
    assert kvp["was_successful"] is True
    assert "error_message" not in kvp

def test_write_and_read_failed_result(temp_db):
    """
    Tests writing a failed DFT result to the database and reading it back.
    """
    db = AseDB(db_path=temp_db)
    atoms = bulk("Ar", "fcc", a=3.8)

    # For a failed calculation, energy/forces/stress might be nonsense or zero
    dft_result = DFTResult(
        total_energy_ev=0.0,
        forces=[[0.0, 0.0, 0.0]] * len(atoms),
        stress=[[0.0] * 3] * 3,
        was_successful=False,
        error_message="SCF failed to converge"
    )

    db_id = db.write(atoms, dft_result)
    retrieved_atoms, kvp = db.get(db_id)

    # Verify key-value pairs for failure case
    assert kvp["was_successful"] is False
    assert kvp["error_message"] == "SCF failed to converge"
