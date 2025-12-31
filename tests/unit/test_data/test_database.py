import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.database import AseDB


@pytest.fixture
def temp_db_path(tmp_path):
    """Fixture to create a temporary database file for each test."""
    db_path = tmp_path / "test.db"
    return str(db_path)


def create_test_atoms():
    """Helper function to create a list of simple Atoms objects."""
    return [
        Atoms("H2", positions=[[0, 0, 0], [0.74, 0, 0]]),
        Atoms("O2", positions=[[0, 0, 0], [1.21, 0, 0]]),
    ]


def test_init_asedb(temp_db_path):
    """Tests that the AseDB class can be initialized."""
    db = AseDB(db_path=temp_db_path)
    assert db.db_path == temp_db_path


def test_init_with_empty_path():
    """Tests that ValueError is raised for an empty path."""
    with pytest.raises(ValueError):
        AseDB(db_path="")


def test_add_structures(temp_db_path):
    """Tests adding structures to the database."""
    db = AseDB(db_path=temp_db_path)
    atoms_list = create_test_atoms()

    db.add_structures(atoms_list, status="test_status")

    # Verify by reading back
    all_data = db.get_all_structures()
    assert len(all_data) == 2
    assert all_data[0][1]["status"] == "test_status"
    assert all_data[1][0].get_chemical_symbols() == ["O", "O"]


def test_get_structures_by_status(temp_db_path):
    """Tests retrieving structures based on their status."""
    db = AseDB(db_path=temp_db_path)
    atoms_list = create_test_atoms()
    db.add_structures(atoms_list[:1], status="status_A")
    db.add_structures(atoms_list[1:], status="status_B")

    atoms_a = db.get_structures_by_status("status_A")
    atoms_b = db.get_structures_by_status("status_B")
    atoms_c = db.get_structures_by_status("status_C")

    assert len(atoms_a) == 1
    assert len(atoms_b) == 1
    assert len(atoms_c) == 0
    assert atoms_a[0].get_chemical_symbols() == ["H", "H"]


def test_update_structure(temp_db_path):
    """Tests updating a structure's status and data."""
    db = AseDB(db_path=temp_db_path)
    atoms_list = create_test_atoms()
    db.add_structures(atoms_list)

    # Get one, "label" it, and update
    h2_atoms = db.get_structures_by_status("needs_labelling")[0]
    assert h2_atoms.get_chemical_symbols() == ["H", "H"]

    # Simulate adding calculator results
    calc = SinglePointCalculator(
        h2_atoms, energy=-30.0, forces=[[0, 0, 0], [0, 0, 0]]
    )
    h2_atoms.calc = calc

    db.update_structure(h2_atoms, status="labelled")

    # Verify the update
    labelled_atoms = db.get_structures_by_status("labelled")
    unlabelled_atoms = db.get_structures_by_status("needs_labelling")

    assert len(labelled_atoms) == 1
    assert len(unlabelled_atoms) == 1
    assert labelled_atoms[0].get_potential_energy() == -30.0


def test_get_all_structures(temp_db_path):
    """Tests retrieving all structures regardless of status."""
    db = AseDB(db_path=temp_db_path)
    atoms_list = create_test_atoms()
    db.add_structures(atoms_list)

    all_structures = db.get_all_structures()
    assert len(all_structures) == 2
