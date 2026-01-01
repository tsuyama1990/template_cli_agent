import ase.db
from ase import Atoms


class AseDBWrapper:
    """A wrapper around the ASE database for managing atomic structures."""

    def __init__(self):
        self.db = None

    def connect(self, path: str):
        """Connects to the ASE database at the given path."""
        self.db = ase.db.connect(path)
        return self.db

    def add_atoms(self, atoms_list: list[Atoms]):
        """Adds a list of Atoms objects to the database."""
        if not self.db:
            raise ConnectionError("Database not connected.")
        for atoms in atoms_list:
            self.db.write(atoms, state="unlabeled")

    def get_unlabeled_rows(self):
        """Retrieves all rows from the database that have the state 'unlabeled'."""
        if not self.db:
            raise ConnectionError("Database not connected.")
        return list(self.db.select(state="unlabeled"))

    def get_labeled_rows(self):
        """Retrieves all rows from the database that have the state 'labeled'."""
        if not self.db:
            raise ConnectionError("Database not connected.")
        return list(self.db.select(state="labeled"))

    def update_row(self, row_id: int, data: dict, key_value_pairs: dict):
        """Updates a specific row in the database with new data and key-value pairs."""
        if not self.db:
            raise ConnectionError("Database not connected.")
        self.db.update(row_id, data=data, key_value_pairs=key_value_pairs)
