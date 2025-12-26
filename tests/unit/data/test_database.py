import pytest
from ase.build import bulk
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult

def test_ase_db_write_and_read(tmp_path):
    """
    Tests writing a successful result to the database and reading it back.
    """
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)

    # 1. Prepare test data
    atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        total_energy_ev=-150.75,
        forces=[[0.0, 0.0, 0.0]] * 2,
        stress=[[0.01] * 3] * 3,
        was_successful=True,
    )

    # 2. Write to the database
    db_id = db.write(atoms, dft_result)
    assert isinstance(db_id, int)
    assert db_id > 0

    # 3. Read it back
    row = db.get(db_id)

    # 4. Verify the data
    assert row.id == db_id
    assert row.natoms == 2
    assert row.key_value_pairs["was_successful"] is True
    assert row.key_value_pairs["total_energy_ev"] == pytest.approx(-150.75)

    # The full Pydantic model is stored in the 'data' field, which ASE
    # automatically decodes from JSON into a dict-like object.
    retrieved_data = row.data
    assert retrieved_data["total_energy_ev"] == pytest.approx(-150.75)
    assert retrieved_data["forces"] == [[0.0, 0.0, 0.0]] * 2


def test_ase_db_write_failed_result(tmp_path):
    """
    Tests writing a failed result to the database.
    """
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)

    atoms = bulk("Au", "fcc", a=4.0)
    dft_result = DFTResult(
        total_energy_ev=0.0,
        forces=[],
        stress=[],
        was_successful=False,
        error_message="Convergence failure",
    )

    db_id = db.write(atoms, dft_result)
    row = db.get(db_id)

    assert row.key_value_pairs["was_successful"] is False
    assert "total_energy_ev" not in row.key_value_pairs

    retrieved_data = row.data
    assert retrieved_data["was_successful"] is False
    assert retrieved_data["error_message"] == "Convergence failure"
