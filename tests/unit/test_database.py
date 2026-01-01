import pytest
import numpy as np
from ase.build import molecule
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.config import DFTResult

def test_add_atoms(temp_db: AseDBWrapper):
    """Test adding an atom and checking its initial state."""
    # The fixture already adds one atom, let's check it
    row = temp_db.get_row(1)
    assert row is not None
    assert row.key_value_pairs["state"] == "unlabeled"

    # Add another one
    co2 = molecule("CO2")
    uid = temp_db.add_atoms(co2, state="custom_state")
    assert uid == 2

    row2 = temp_db.get_row(2)
    assert row2.key_value_pairs["state"] == "custom_state"

def test_get_atoms_by_id(temp_db: AseDBWrapper):
    """Test retrieving an atom by its ID."""
    atoms = temp_db.get_atoms_by_id(1)
    assert atoms is not None
    assert len(atoms) == 3 # H2O molecule
    assert sorted(atoms.get_chemical_symbols()) == sorted(['H', 'H', 'O'])

def test_get_nonexistent_atom(temp_db: AseDBWrapper):
    """Test that retrieving a non-existent ID returns None."""
    atoms = temp_db.get_atoms_by_id(999)
    assert atoms is None

def test_update_labels(temp_db: AseDBWrapper):
    """Test updating a row with DFT results."""
    dft_result = DFTResult(
        energy=-450.0,
        forces=np.array([[0.1, 0.0, 0.0]] * 3),
        stress=np.eye(3) * 0.01
    )

    temp_db.update_labels(1, dft_result)

    row = temp_db.get_row(1)
    assert row.key_value_pairs["state"] == "labeled"

    # Verify the stored JSON data
    stored_json = row.key_value_pairs["dft_result"]
    reloaded_result = DFTResult.model_validate_json(stored_json)

    assert reloaded_result.energy == dft_result.energy
    np.testing.assert_array_equal(reloaded_result.forces, dft_result.forces)

def test_get_all_labeled_atoms(temp_db: AseDBWrapper):
    """Test retrieving all labeled atoms from the database."""
    # Initially, no labeled atoms
    labeled = temp_db.get_all_labeled_atoms()
    assert len(labeled) == 0

    # Label the first atom
    dft_result1 = DFTResult(energy=-1.0, forces=np.zeros((3, 3)), stress=np.zeros((3, 3)))
    temp_db.update_labels(1, dft_result1)

    # Add and label a second atom
    co2 = molecule("CO2")
    uid2 = temp_db.add_atoms(co2)
    dft_result2 = DFTResult(energy=-2.0, forces=np.ones((3, 3)), stress=np.ones((3, 3)))
    temp_db.update_labels(uid2, dft_result2)

    labeled = temp_db.get_all_labeled_atoms()
    assert len(labeled) == 2

    # Check the contents
    symbols = [sorted(res[0].get_chemical_symbols()) for res in labeled]
    energies = [res[1].energy for res in labeled]

    assert sorted(['H', 'H', 'O']) in symbols
    assert sorted(['C', 'O', 'O']) in symbols
    assert -1.0 in energies
    assert -2.0 in energies
