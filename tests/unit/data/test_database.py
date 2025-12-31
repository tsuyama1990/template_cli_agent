import os

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.database import AseDBWrapper


@pytest.fixture
def temp_db():
    db_path = "test.db"
    yield db_path
    if os.path.exists(db_path):
        os.remove(db_path)


def test_add_and_get_atoms(temp_db):
    db = AseDBWrapper(temp_db)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    id = db.add_atoms(atoms, labeled=False)

    retrieved_row = db.get_atoms_by_id(id)
    assert retrieved_row.id == id
    assert len(retrieved_row.toatoms()) == 2


def test_get_rows_to_label(temp_db):
    db = AseDBWrapper(temp_db)
    atoms1 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    atoms2 = Atoms("O2", positions=[[0, 0, 0], [0, 0, 1.2]])

    db.add_atoms(atoms1, labeled=False)

    id2 = db.add_atoms(atoms2, labeled=False)
    calc = SinglePointCalculator(atoms2, energy=-1.0, forces=[[0, 0, 0], [0, 0, 0]])
    atoms2.calc = calc
    db.update_row_with_dft_results(id2, atoms2)

    rows_to_label = db.get_rows_to_label()
    assert len(rows_to_label) == 1
    assert rows_to_label[0].formula == "H2"


def test_update_row_with_dft_results(temp_db):
    db = AseDBWrapper(temp_db)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    id = db.add_atoms(atoms)

    calc = SinglePointCalculator(atoms, energy=-1.0, forces=[[0, 0, 0], [0, 0, 0]])
    atoms.calc = calc
    db.update_row_with_dft_results(id, atoms)

    updated_row = db.get_atoms_by_id(id)
    assert updated_row.labeled is True
    assert updated_row.energy == -1.0
