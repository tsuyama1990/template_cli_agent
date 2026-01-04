
from ase import Atoms
from ase.db import connect


class DatabaseManager:
    def __init__(self, db_path: str):
        self.db_path = db_path

    def write_atoms(self, atoms_list: list[Atoms]):
        """Writes a list of Atoms objects to the database."""
        with connect(self.db_path) as db:
            for atoms in atoms_list:
                db.write(atoms)
