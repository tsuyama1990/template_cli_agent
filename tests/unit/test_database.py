# tests/unit/test_database.py

from unittest.mock import MagicMock, call

import pytest
from ase import Atoms

from mlip_autopipec.data.database import AseDBWrapper


@pytest.fixture
def mock_db_connection():
    """Fixture to create a mock ASE database connection."""
    # Create a mock for the row object
    mock_row = MagicMock()
    mock_row.toatoms.return_value = Atoms("H")
    mock_row.key_value_pairs = {"state": "labeled"}
    mock_row.data = {"energy": -1.0, "forces": [[0, 0, 0]]}  # Add mock data

    # Create a mock for the connection object
    mock_conn = MagicMock()
    mock_conn.select.return_value = [mock_row]
    mock_conn.write.return_value = 1
    mock_conn.update.return_value = None
    return mock_conn


@pytest.fixture
def mock_ase_db_connect(mocker, mock_db_connection):
    """Fixture to mock the ase.db.connect context manager."""
    mock_connect = mocker.patch("ase.db.connect")
    # The __enter__ method of the context manager should return the mock connection
    mock_connect.return_value.__enter__.return_value = mock_db_connection
    return mock_connect


def test_add_atoms(mock_ase_db_connect, mock_db_connection):
    """Tests that atoms are added to the database with the correct state."""
    db_wrapper = AseDBWrapper("test.db")
    atoms1 = Atoms("H")
    atoms2 = Atoms("He")

    db_wrapper.add_atoms([atoms1, atoms2])

    mock_ase_db_connect.assert_called_once_with("test.db")
    # Check that write was called for each atom
    assert mock_db_connection.write.call_count == 2
    mock_db_connection.write.assert_has_calls(
        [call(atoms1, state="unlabeled"), call(atoms2, state="unlabeled")]
    )


def test_get_unlabeled_rows(mock_ase_db_connect, mock_db_connection):
    """Tests that unlabeled rows are correctly queried from the database."""
    db_wrapper = AseDBWrapper("test.db")

    rows = db_wrapper.get_unlabeled_rows()

    mock_ase_db_connect.assert_called_once_with("test.db")
    mock_db_connection.select.assert_called_once_with(state="unlabeled")
    assert len(rows) == 1  # As per mock_db_connection setup
    assert rows[0] is mock_db_connection.select.return_value[0]


def test_update_row(mock_ase_db_connect, mock_db_connection):
    """Tests that a database row is updated with the correct data."""
    db_wrapper = AseDBWrapper("test.db")
    row_id = 1
    data = {"energy": -1.0, "forces": [[0, 0, 0]]}
    kvp = {"state": "labeled"}

    db_wrapper.update_row(row_id, data, kvp)

    mock_ase_db_connect.assert_called_once_with("test.db")
    mock_db_connection.update.assert_called_once_with(
        row_id, data=data, state="labeled"
    )


def test_get_labeled_atoms(mock_ase_db_connect, mock_db_connection):
    """Tests that labeled atoms are correctly retrieved and parsed."""
    db_wrapper = AseDBWrapper("test.db")

    labeled_data = db_wrapper.get_labeled_atoms()

    mock_ase_db_connect.assert_called_once_with("test.db")
    mock_db_connection.select.assert_called_once_with(state="labeled")
    assert len(labeled_data) == 1
    atoms, kvp = labeled_data[0]
    assert isinstance(atoms, Atoms)
    # Corrected assertion to check for the combined dictionary
    expected_kvp = {
        "state": "labeled",
        "data": {"energy": -1.0, "forces": [[0, 0, 0]]},
    }
    assert kvp == expected_kvp
    mock_db_connection.select.return_value[0].toatoms.assert_called_once()
