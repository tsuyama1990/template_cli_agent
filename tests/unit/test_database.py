import os
import pytest
from ase.db import connect
from ase.build import bulk
from mlip_autopipec.data.database import AseDB
from ase.calculators.singlepoint import SinglePointCalculator

@pytest.fixture
def db_path():
    db_file = "test_db.db"
    yield db_file
    if os.path.exists(db_file):
        os.remove(db_file)

def test_asedb_write_and_read(db_path):
    db = AseDB(db_path)
    atoms = bulk("Si", "diamond", a=5.43)

    # Attach data via a calculator
    calc = SinglePointCalculator(atoms, energy=-100.0, forces=[[0.0, 0.0, 0.0]] * len(atoms))
    atoms.calc = calc

    # Write with metadata
    metadata = {'was_successful': True}
    db_id = db.write(atoms, metadata=metadata)

    assert isinstance(db_id, int)

    with connect(db_path) as con:
        assert len(con) == 1
        row = con.get(id=db_id)
        assert row.energy == -100.0
        assert row.key_value_pairs['was_successful'] is True

def test_asedb_invalid_path():
    with pytest.raises(ValueError):
        AseDB("test_db.txt")
