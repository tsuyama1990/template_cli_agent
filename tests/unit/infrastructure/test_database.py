"""Unit tests for the AseDB infrastructure adapter."""

import pytest
from ase import Atoms
from ase.db import connect
import numpy as np
import os

from mlip_autopipec.infrastructure.database import AseDB
from mlip_autopipec.domain.models import DFTResult

@pytest.fixture
def db_path(tmp_path):
    """Pytest fixture to create a temporary database for testing."""
    path = os.path.join(tmp_path, "test.db")
    yield path
    if os.path.exists(path):
        os.remove(path)

def test_add_atoms(db_path):
    db = AseDB(db_path)
    atoms1 = Atoms('H')
    db.add_atoms([atoms1])
    conn = connect(db_path)
    assert len(conn) == 1
    assert conn.get(id=1).state == 'unlabelled'

def test_get_atoms_by_state(db_path):
    db = AseDB(db_path)
    db.add_atoms([Atoms('H')])
    conn = connect(db_path)
    conn.write(Atoms('O'), state='labelled')
    unlabelled = db.get_atoms_by_state('unlabelled')
    labelled = db.get_atoms_by_state('labelled')
    assert len(unlabelled) == 1
    assert unlabelled[0].get_chemical_symbols() == ['H']
    assert len(labelled) == 1

def test_update_with_dft_results(db_path):
    db = AseDB(db_path)
    db.add_atoms([Atoms('He', cell=np.eye(3), pbc=True)])
    conn = connect(db_path)
    atoms_id = conn.get(state='unlabelled').id
    stress_3x3 = np.array([[0.1, 0.4, 0.5], [0.4, 0.2, 0.6], [0.5, 0.6, 0.3]])
    dft_results = DFTResult(
        energy=-1.0,
        forces=np.array([[0.1, 0.2, 0.3]]),
        stress=stress_3x3
    )
    db.update_with_dft_results(atoms_id, dft_results)
    updated_row = conn.get(id=atoms_id)
    updated_atoms = updated_row.toatoms()
    assert updated_row.state == 'labelled'
    assert updated_atoms.get_potential_energy() == -1.0
    np.testing.assert_array_almost_equal(updated_atoms.get_stress(voigt=False), dft_results.stress)
