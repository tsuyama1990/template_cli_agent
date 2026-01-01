from unittest.mock import MagicMock, patch, call

import pytest
from ase import Atoms

# The class we will implement in Phase 3
# We are writing tests for it now (TDD)
from mlip_autopipec.data.database import AseDBWrapper


@pytest.fixture
def db_wrapper():
    """Provides an instance of AseDBWrapper."""
    return AseDBWrapper()


# Patch where 'connect' is LOOKED UP, which is in the database module
@patch("mlip_autopipec.data.database.ase.db.connect")
def test_connect(mock_ase_connect, db_wrapper):
    """Test that the connect method calls ase.db.connect and sets the db attribute."""
    mock_connection = MagicMock()
    mock_ase_connect.return_value = mock_connection

    db_wrapper.connect("test.db")

    mock_ase_connect.assert_called_once_with("test.db")
    assert db_wrapper.db is mock_connection


def test_add_atoms(db_wrapper):
    """Test adding a list of Atoms objects."""
    db_wrapper.db = MagicMock()  # Simulate a connection
    atoms1 = Atoms("H", positions=[[0, 0, 0]])
    atoms2 = Atoms("He", positions=[[1, 1, 1]])

    db_wrapper.add_atoms([atoms1, atoms2])

    assert db_wrapper.db.write.call_count == 2
    # Use any_call since the order doesn't matter
    # CORRECTED: Added state='unlabeled' to the assertion
    db_wrapper.db.write.assert_any_call(atoms1, state="unlabeled")
    db_wrapper.db.write.assert_any_call(atoms2, state="unlabeled")


def test_add_atoms_no_connection_raises_error(db_wrapper):
    """Test that methods fail if the database is not connected."""
    with pytest.raises(ConnectionError, match="Database not connected"):
        db_wrapper.add_atoms([])


def test_get_unlabeled_rows(db_wrapper):
    """Test retrieving rows that are marked as 'unlabeled'."""
    db_wrapper.db = MagicMock()
    mock_row = MagicMock()
    mock_row.toatoms.return_value = Atoms("Li")
    db_wrapper.db.select.return_value = [mock_row]

    rows = db_wrapper.get_unlabeled_rows()

    db_wrapper.db.select.assert_called_once_with(state="unlabeled")
    assert len(rows) == 1
    assert rows[0] is mock_row


def test_update_row(db_wrapper):
    """Test updating a row with new data and key-value pairs."""
    db_wrapper.db = MagicMock()
    row_id = 1
    data = {"energy": -123.45}
    kvp = {"state": "labeled", "calculation_time": 100.0}

    db_wrapper.update_row(row_id, data=data, key_value_pairs=kvp)

    db_wrapper.db.update.assert_called_once_with(
        row_id, data=data, key_value_pairs=kvp
    )
