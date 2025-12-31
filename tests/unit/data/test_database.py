
import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from mlip_autopipec.data.database import AseDBWrapper


@pytest.fixture
def temp_db_path(tmp_path):
    """Fixture to create a temporary database path."""
    db_file = tmp_path / "test.db"
    yield str(db_file)
    if db_file.exists():
        db_file.unlink()

@pytest.fixture
def db_wrapper(temp_db_path):
    """Fixture to create an AseDBWrapper instance."""
    return AseDBWrapper(temp_db_path)

def test_init_raises_error_for_empty_path():
    """Test that ValueError is raised if db_path is empty."""
    with pytest.raises(ValueError):
        AseDBWrapper("")

def test_add_atoms(db_wrapper, temp_db_path):
    """Test adding atoms to the database."""
    atoms1 = Atoms('H', positions=[(0, 0, 0)])
    atoms2 = Atoms('He', positions=[(1, 1, 1)])

    db_wrapper.add_atoms([atoms1, atoms2], source="test")

    with connect(temp_db_path) as db:
        assert len(db) == 2
        row1 = db.get(1)
        row2 = db.get(2)
        assert "".join(row1.symbols) == 'H'
        assert "".join(row2.symbols) == 'He'
        assert row1.key_value_pairs['labelled'] is False
        assert row2.key_value_pairs['labelled'] is False
        assert row1.key_value_pairs['source'] == "test"
        assert row2.key_value_pairs['source'] == "test"

def test_get_row(db_wrapper):
    """Test retrieving a single row by ID."""
    atoms = Atoms('C', positions=[(0, 0, 0)])
    db_wrapper.add_atoms([atoms])

    row = db_wrapper.get_row(1)
    assert row is not None
    assert "".join(row.symbols) == 'C'

    non_existent_row = db_wrapper.get_row(999)
    assert non_existent_row is None

def test_get_rows_to_label(db_wrapper, temp_db_path):
    """Test retrieving only unlabeled rows."""
    atoms_unlabeled = Atoms('N')
    atoms_labeled = Atoms('O')

    db_wrapper.add_atoms([atoms_unlabeled]) # Adds with labelled=False
    with connect(temp_db_path) as db:
        db.write(atoms_labeled, key_value_pairs={'labelled': True})

    rows_to_label = db_wrapper.get_rows_to_label()
    assert len(rows_to_label) == 1
    assert "".join(rows_to_label[0].symbols) == 'N'
    assert rows_to_label[0].key_value_pairs['labelled'] is False

def test_update_row_with_dft_results(db_wrapper):
    """Test updating a row with DFT calculation results."""
    atoms = Atoms('Si', positions=[(0, 0, 0)])
    db_wrapper.add_atoms([atoms])

    dft_results = {
        'energy': -10.0,
        'forces': np.array([[0.1, 0.2, 0.3]]),
        'stress': np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]), # Voigt form
    }

    db_wrapper.update_row_with_dft_results(1, dft_results)

    updated_row = db_wrapper.get_row(1)
    assert updated_row is not None
    assert updated_row.key_value_pairs['labelled'] is True

    # Check that calculator results are correctly stored
    updated_atoms = updated_row.toatoms()
    assert updated_atoms.get_potential_energy() == pytest.approx(-10.0)
    assert np.allclose(updated_atoms.get_forces(), dft_results['forces'])
    assert np.allclose(updated_atoms.get_stress(), dft_results['stress'])

def test_get_all_labeled_rows(db_wrapper, temp_db_path):
    """Test retrieving all successfully labeled rows."""
    atoms_unlabeled = Atoms('F')
    atoms_labeled = Atoms('Ne')

    db_wrapper.add_atoms([atoms_unlabeled])
    with connect(temp_db_path) as db:
        db.write(atoms_labeled, key_value_pairs={'labelled': True})

    labeled_rows = db_wrapper.get_all_labeled_rows()
    assert len(labeled_rows) == 1
    assert "".join(labeled_rows[0].symbols) == 'Ne'
    assert labeled_rows[0].key_value_pairs['labelled'] is True
