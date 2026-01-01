from typing import List, Tuple, Optional
import ase.db
from ase.atoms import Atoms
from .config import DFTResult

class AseDBWrapper:
    """A wrapper class for the ASE database to manage data persistence."""

    def __init__(self, db_path: str):
        """
        Initializes the database wrapper.

        Args:
            db_path: Path to the SQLite database file.
        """
        self.db_path = db_path

    def _connect(self):
        """Creates a connection to the database."""
        return ase.db.connect(self.db_path)

    def add_atoms(self, atoms: Atoms, state: str = "unlabeled") -> int:
        """
        Adds an Atoms object to the database.

        Args:
            atoms: The ASE Atoms object to add.
            state: The initial state of the structure.

        Returns:
            The ID of the newly inserted row.
        """
        with self._connect() as con:
            uid = con.write(atoms, key_value_pairs={"state": state})
        return uid

    def get_atoms_by_id(self, uid: int) -> Optional[Atoms]:
        """
        Retrieves a single Atoms object by its unique ID.

        Args:
            uid: The unique ID of the row.

        Returns:
            The Atoms object, or None if not found.
        """
        with self._connect() as con:
            try:
                row = con.get(id=uid)
                return row.toatoms()
            except KeyError:
                return None

    def update_labels(self, uid: int, dft_result: DFTResult):
        """
        Updates a database entry with DFT calculation results.

        Args:
            uid: The unique ID of the row to update.
            dft_result: The Pydantic DFTResult object containing the labels.
        """
        with self._connect() as con:
            # ASE's update method uses key=value for kvp updates
            con.update(
                id=uid,
                state="labeled",
                dft_result=dft_result.model_dump_json()
            )

    def get_all_labeled_atoms(self) -> List[Tuple[Atoms, DFTResult]]:
        """
        Retrieves all atoms that are marked as 'labeled'.

        Returns:
            A list of tuples, where each tuple contains the Atoms object
            and the corresponding deserialized DFTResult object.
        """
        results = []
        with self._connect() as con:
            for row in con.select(state="labeled"):
                atoms = row.toatoms()
                dft_result_json = row.key_value_pairs["dft_result"]
                dft_result = DFTResult.model_validate_json(dft_result_json)
                results.append((atoms, dft_result))
        return results

    def get_row(self, uid: int):
        """Gets a raw database row by ID."""
        with self._connect() as con:
            try:
                return con.get(id=uid)
            except KeyError:
                return None
