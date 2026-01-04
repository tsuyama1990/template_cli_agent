from ase import Atoms
from ase.db import connect


class DatabaseManager:
    """
    Handles all interactions with the ASE database.
    """

    def __init__(self, db_path: str) -> None:
        self.db_path = db_path

    def write_structures(self, structures: list[Atoms]):
        """
        Writes a list of Atoms objects to the database.
        """
        with connect(self.db_path) as db:
            for atoms in structures:
                db.write(atoms)
