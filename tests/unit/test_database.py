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
    dft_result = DFTResult(
        energy=-50.0, forces=np.zeros((1, 3)), stress=np.zeros((3, 3))
    )
    db_wrapper.update_labels(id_labeled, dft_result)

    labeled_atoms = db_wrapper.get_all_labeled_atoms()
    assert len(labeled_atoms) == 1
    assert labeled_atoms[0][0].get_chemical_symbols() == ["Be"]
