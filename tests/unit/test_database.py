import pytest
from ase.build import bulk
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult

def test_db_write_and_read(tmp_path):
    """
    Tests that data can be written to and read from the AseDB.
    """
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)

    # 1. Create sample data
    atoms = bulk("Si", "diamond", a=5.43)
    dft_result = DFTResult(
        total_energy_ev=-150.0,
        forces=[[0.0, 0.0, 0.0]] * 2,
        stress=[[0.0] * 3] * 2,
        was_successful=True,
    )

    # 2. Write data to the database
    db_id = db.write(atoms, dft_result)
    assert isinstance(db_id, int)
    assert db_id > 0

    # 3. Read data back from the database
    retrieved_data = db.get_row(db_id)
    retrieved_atoms = db.get_atoms(db_id)

    # 4. Verify the retrieved data
    assert retrieved_data is not None
    assert retrieved_atoms is not None

    # Check that the DFTResult data matches
    retrieved_dft_result = DFTResult(**retrieved_data)
    assert retrieved_dft_result.was_successful is True
    assert retrieved_dft_result.total_energy_ev == pytest.approx(-150.0)
    assert len(retrieved_dft_result.forces) == 2

    # Check that the Atoms object is reconstructed correctly
    assert len(retrieved_atoms) == 2
    assert retrieved_atoms.get_chemical_symbols() == ["Si", "Si"]
    assert retrieved_atoms.cell.bandpath().path == 'GXWKGLUWLK,UX' # Verifies it's a diamond lattice

def test_get_non_existent_row(tmp_path):
    """
    Tests that getting a non-existent row returns None.
    """
    db_path = tmp_path / "test.db"
    db = AseDB(db_path)

    assert db.get_row(999) is None
    assert db.get_atoms(999) is None
