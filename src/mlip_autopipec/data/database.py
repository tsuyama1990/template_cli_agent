
import ase.db
from ase import Atoms


class AseDB:
    """A wrapper class for the ASE database to provide a domain-specific API."""

    def __init__(self, db_path: str):
        """
        Initializes the database connection.

        Args:
            db_path: Path to the ASE database file.
        """
        if not db_path:
            raise ValueError("Database path cannot be empty.")
        self.db_path = db_path

    def _connect(self):
        """Creates a connection to the database."""
        return ase.db.connect(self.db_path)

    def add_structures(self, structures: list[Atoms], status: str = "needs_labelling"):
        """
        Adds a list of Atoms objects to the database.

        Args:
            structures: A list of ase.Atoms objects to add.
            status: The initial status to assign to the structures.
        """
        with self._connect() as db:
            for atoms in structures:
                db.write(atoms, key_value_pairs={"status": status})

    def get_structures_by_status(self, status: str) -> list[Atoms]:
        """
        Retrieves all structures with a given status.

        Args:
            status: The status to filter by (e.g., 'needs_labelling').

        Returns:
            A list of matching ase.Atoms objects.
        """
        atoms_list = []
        with self._connect() as db:
            for row in db.select(f"status={status}"):
                atoms_list.append(row.toatoms())
        return atoms_list

    def update_structure(self, atoms: Atoms, status: str):
        """
        Updates an existing structure in the database, identified by its unique ID.
        This is a bit tricky with ase.db, as direct updates are based on row id.
        A more robust way is to find the row by a unique property if available,
        or just re-write it if the properties change. For now, we assume we
        can select it and update its metadata.

        Args:
            atoms: The ase.Atoms object with updated information (e.g., calculator).
            status: The new status to set.
        """
        # ASE DB doesn't have a great way to find a specific atoms object
        # without a unique ID. We'll assume for now we can select based on
        # a unique property if we had one. A simpler approach for this project
        # is to get the row id and then update.
        # This implementation is a placeholder for a more robust one if needed.
        with self._connect() as db:
            # This is not efficient, but ase.db lacks a direct update-by-content.
            # We find the row ID of the structure that is geometrically identical.
            # This is fragile and assumes no two identical structures exist.
            # In a real scenario, we'd manage unique IDs.
            found_id = None
            for row in db.select():
                if (
                    len(atoms) == len(row.toatoms())
                    and (atoms.get_positions() == row.toatoms().get_positions()).all()
                ):
                    found_id = row.id
                    break

            if found_id:
                db.update(found_id, atoms=atoms, status=status)
            else:
                # If not found, perhaps it's a new one. This logic might need refinement.
                self.add_structures([atoms], status=status)

    def get_all_structures(self) -> list[tuple[Atoms, dict]]:
        """Gets all structures and their key-value pairs."""
        results = []
        with self._connect() as db:
            for row in db.select():
                results.append((row.toatoms(), row.key_value_pairs))
        return results
