from pathlib import Path

from ase import Atoms
from ase.db import connect

from mlip_autopipec.config import DFTResult


class AseDBWrapper:
    """
    A wrapper class for the ASE database to manage atomic structures and their labels.

    This class provides a high-level API for interacting with the ASE database,
    handling the serialization and deserialization of DFT results and managing
    the state of each structure in the database.

    Note on SQL Injection: This class relies on the `ase.db` library, which uses
    parameterized queries under the hood, thus protecting against SQL injection.
    """

    def __init__(self, db_path: str):
        """
        Initializes the database wrapper.

        Args:
            db_path: Path to the ASE database file.
        """
        self.db_path = Path(db_path)

    def _connect(self):
        """Connects to the ASE database."""
        return connect(self.db_path)

    def add_atoms(
        self, atoms: Atoms, state: str = "unlabeled", is_initial: bool = False
    ) -> int:
        """Adds an Atoms object to the database.

        Args:
            atoms: The ASE Atoms object to add.
            state: The initial state of the atoms object.
            is_initial: Flag to mark if this is an initial structure.

        Returns:
            The ID of the newly inserted row.
        """
        with self._connect() as db:
            kvp = {"state": state, "is_initial": is_initial}
            id = db.write(atoms, key_value_pairs=kvp)
        return id

    def get_all_initial_structures(self) -> list[Atoms]:
        """Retrieves all structures marked as 'initial'.

        Returns:
            A list of Atoms objects.
        """
        initial_structures = []
        with self._connect() as db:
            # The select method in ASE db treats boolean False as not present
            # so we select for `is_initial=True`
            for row in db.select(is_initial=True):
                initial_structures.append(row.toatoms())
        return initial_structures

    def get_atoms_by_id(self, id: int) -> Atoms:
        """Retrieves an Atoms object by its ID.

        Args:
            id: The ID of the Atoms object to retrieve.

        Returns:
            The retrieved Atoms object.
        """
        with self._connect() as db:
            row = db.get(id=id)
            return row.toatoms()

    def update_labels(self, id: int, dft_result: DFTResult):
        """Updates an entry with DFT calculation results.

        Args:
            id: The ID of the entry to update.
            dft_result: The DFTResult object containing the calculated labels.
        """
        with self._connect() as db:
            dft_result_json = dft_result.model_dump_json()
            db.update(
                id,
                state="labeled",
                dft_result=dft_result_json,
            )

    def update_state(self, id: int, state: str):
        """Updates the state of an entry.

        Args:
            id: The ID of the entry to update.
            state: The new state to set.
        """
        with self._connect() as db:
            db.update(id, state=state)

    def get_all_labeled_atoms(self) -> list[tuple[Atoms, DFTResult]]:
        """Retrieves all atoms that have been successfully labeled.

        Returns:
            A list of tuples, each containing an Atoms object and its corresponding DFTResult.
        """
        results = []
        with self._connect() as db:
            for row in db.select(state="labeled"):
                try:
                    atoms = row.toatoms()
                    dft_result_json = row.key_value_pairs["dft_result"]
                    dft_result = DFTResult.model_validate_json(dft_result_json)
                    results.append((atoms, dft_result))
                except Exception as e:
                    import logging

                    logging.warning(f"Skipping corrupted data in database (ID: {row.id}): {e}")
        return results

    def get_unlabeled_ids(self) -> list[int]:
        """Retrieves the IDs of all unlabeled atoms.

        Returns:
            A list of integer IDs.
        """
        ids = []
        with self._connect() as db:
            for row in db.select(state="unlabeled"):
                ids.append(row.id)
        return ids

    def is_empty(self) -> bool:
        """Checks if the database is empty.

        Returns:
            True if the database has no entries, False otherwise.
        """
        with self._connect() as db:
            return db.count() == 0
