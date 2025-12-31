"""Unit tests for the AseDB wrapper."""

import os

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import DFTResult


@pytest.fixture
def db_path(tmp_path):
    """Pytest fixture to create a temporary database for testing."""
    path = os.path.join(tmp_path, "test.db")
    yield path
    if os.path.exists(path):
        os.remove(path)


def test_add_atoms(db_path):
    """Test adding atoms to the database."""
    db = AseDB(db_path)
    atoms1 = Atoms("H", positions=[(0, 0, 0)])
    atoms2 = Atoms("O", positions=[(0, 0, 0)])
    db.add_atoms([atoms1, atoms2])

    # Verify directly with ase.db.connect
    conn = connect(db_path)
    assert len(conn) == 2
    assert conn.get(id=1).state == "unlabelled"
    assert conn.get(id=2).state == "unlabelled"


def test_get_atoms_by_state(db_path):
    """Test retrieving atoms by their state."""
    db = AseDB(db_path)
    atoms_unlabelled = [Atoms("H")]
    atoms_labelled = [Atoms("O")]

    db.add_atoms(atoms_unlabelled)
    # Manually update one to 'labelled' for the test
    conn = connect(db_path)
    conn.write(atoms_labelled[0], state="labelled")

    unlabelled = db.get_atoms_by_state("unlabelled")
    labelled = db.get_atoms_by_state("labelled")

    assert len(unlabelled) == 1
    assert unlabelled[0].get_chemical_symbols() == ["H"]
    assert unlabelled[0].info["db_id"] is not None  # check id is attached
    assert len(labelled) == 1
    assert labelled[0].get_chemical_symbols() == ["O"]


def test_update_with_dft_results(db_path):
    """Test updating an entry with DFT results."""
    db = AseDB(db_path)
    atoms = Atoms("He", positions=[(0, 0, 0)], cell=np.eye(3), pbc=True)
    db.add_atoms([atoms])

    # Get the id of the added atom
    conn = connect(db_path)
    row = conn.get(state="unlabelled")
    atoms_id = row.id

    stress_3x3 = np.array([[0.1, 0.4, 0.5], [0.4, 0.2, 0.6], [0.5, 0.6, 0.3]])

    dft_results = DFTResult(energy=-1.0, forces=np.array([[0.1, 0.2, 0.3]]), stress=stress_3x3)

    db.update_with_dft_results(atoms_id, dft_results)

    updated_row = conn.get(id=atoms_id)
    updated_atoms = updated_row.toatoms()

    assert updated_row.state == "labelled"
    assert updated_atoms.get_potential_energy() == -1.0
    np.testing.assert_array_almost_equal(updated_atoms.get_forces(), dft_results.forces)
    # ase.db stores stress in voigt form. get_stress(voigt=False) converts it back.
    np.testing.assert_array_almost_equal(updated_atoms.get_stress(voigt=False), dft_results.stress)
