import tempfile
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.config import DFTResult
from mlip_autopipec.database import AseDBWrapper


@pytest.fixture
def temp_db():
    """Creates a temporary database for testing."""
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = f.name
    yield db_path
    Path(db_path).unlink()


def test_add_atoms(temp_db):
    """Tests adding an Atoms object to the database."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    atoms = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    id = db_wrapper.add_atoms(atoms)
    assert id == 1

    retrieved_atoms = db_wrapper.get_atoms_by_id(id)
    assert len(retrieved_atoms) == 3
    assert retrieved_atoms.get_chemical_symbols() == ["H", "H", "O"]


def test_update_labels(temp_db):
    """Tests updating an entry with DFT labels."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    atoms = Atoms("Si", positions=[[0, 0, 0]])
    id = db_wrapper.add_atoms(atoms)

    dft_result = DFTResult(
        energy=-100.0,
        forces=np.ones((1, 3)),
        stress=np.eye(3) * 2,
    )

    db_wrapper.update_labels(id, dft_result)

    labeled_data = db_wrapper.get_all_labeled_atoms()
    assert len(labeled_data) == 1

    retrieved_atoms, retrieved_result = labeled_data[0]
    assert retrieved_atoms.get_chemical_symbols() == ["Si"]
    assert retrieved_result.energy == dft_result.energy
    np.testing.assert_array_almost_equal(retrieved_result.forces, dft_result.forces)
    np.testing.assert_array_almost_equal(retrieved_result.stress, dft_result.stress)


def test_get_all_labeled_atoms(temp_db):
    """Tests retrieving all labeled atoms from the database."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    db_wrapper.add_atoms(Atoms("H"), state="unlabeled")
    db_wrapper.add_atoms(Atoms("Li"), state="labeling_failed")

    # Manually add a labeled atom to check retrieval
    id_labeled = db_wrapper.add_atoms(Atoms("Be"))
    dft_result = DFTResult(energy=-50.0, forces=np.zeros((1, 3)), stress=np.zeros((3, 3)))
    db_wrapper.update_labels(id_labeled, dft_result)

    labeled_atoms = db_wrapper.get_all_labeled_atoms()
    assert len(labeled_atoms) == 1
    assert labeled_atoms[0][0].get_chemical_symbols() == ["Be"]


def test_get_all_labeled_atoms_corrupted_data(temp_db):
    """Tests that corrupted data in the database is gracefully skipped."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    atoms = Atoms("C")
    id_corrupted = db_wrapper.add_atoms(atoms)

    # Manually insert corrupted data
    with db_wrapper._connect() as db:
        db.update(id_corrupted, state="labeled", dft_result="invalid json")

    labeled_atoms = db_wrapper.get_all_labeled_atoms()
    assert len(labeled_atoms) == 0


def test_database_connection_error(tmp_path):
    """Tests that a connection error is handled gracefully."""
    from sqlite3 import OperationalError

    db_path = tmp_path / "non_existent_dir" / "test.db"
    db_wrapper = AseDBWrapper(db_path=str(db_path))
    with pytest.raises(OperationalError):
        db_wrapper.add_atoms(Atoms("H"))


def test_get_unlabeled_ids(temp_db):
    """Tests retrieving the IDs of all unlabeled atoms."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    id1 = db_wrapper.add_atoms(Atoms("H"), state="unlabeled")
    id2 = db_wrapper.add_atoms(Atoms("He"), state="unlabeled")
    db_wrapper.add_atoms(Atoms("Li"), state="labeled")

    unlabeled_ids = db_wrapper.get_unlabeled_ids()
    assert set(unlabeled_ids) == {id1, id2}


def test_is_empty(temp_db):
    """Tests the is_empty method."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    assert db_wrapper.is_empty()

    db_wrapper.add_atoms(Atoms("H"))
    assert not db_wrapper.is_empty()


def test_update_labels_sets_state_correctly(temp_db):
    """Tests that update_labels correctly sets the state and DFT result."""
    db_wrapper = AseDBWrapper(db_path=temp_db)
    atoms = Atoms("Si", positions=[[0, 0, 0]])
    id = db_wrapper.add_atoms(atoms, state="unlabeled")

    dft_result = DFTResult(
        energy=-100.0,
        forces=np.ones((1, 3)),
        stress=np.eye(3) * 2,
    )

    db_wrapper.update_labels(id, dft_result)

    # Directly query the database to verify the state and data
    with db_wrapper._connect() as db:
        row = db.get(id=id)
        assert row.key_value_pairs["state"] == "labeled"
        retrieved_result = DFTResult.model_validate_json(row.key_value_pairs["dft_result"])
        assert retrieved_result.energy == dft_result.energy
