from ase.atoms import Atoms
from ase.db import connect
from ase.db.row import AtomsRow


class AseDBWrapper:
    """A wrapper for the ASE database to handle data persistence."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.connection = connect(self.db_path)

    def add_atoms(self, atoms: Atoms, **kwargs) -> int:
        """Adds a new Atoms object to the database."""
        return self.connection.write(atoms, **kwargs)

    def get_atoms_by_id(self, id: int) -> AtomsRow:
        """Retrieves a single Atoms object by its ID."""
        return self.connection.get(id=id)

    def get_rows_to_label(self) -> list[AtomsRow]:
        """Retrieves all rows that have not yet been labeled."""
        return list(self.connection.select(labeled=False))

    def update_row_with_dft_results(self, id: int, atoms: Atoms):
        """Updates a row with a new Atoms object containing calculator results."""
        self.connection.update(id, atoms=atoms, labeled=True)
