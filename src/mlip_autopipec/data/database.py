# src/mlip_autopipec/data/database.py

from typing import Any

import ase.db
from ase import Atoms


class AseDBWrapper:
    """A wrapper around the ASE database for managing atomic structures."""

    def __init__(self, db_path: str):
        """
        Initializes the database wrapper.

        Args:
            db_path: The path to the ASE database file.
        """
        self.db_path = db_path

    def add_atoms(self, atoms_list: list[Atoms]):
        """
        Adds a list of ASE Atoms objects to the database.

        Args:
            atoms_list: A list of Atoms objects to add.
        """
        with ase.db.connect(self.db_path) as db:
            for atoms in atoms_list:
                db.write(atoms, state="unlabeled")

    def get_unlabeled_rows(self) -> list[Any]:
        """
        Retrieves all rows from the database that are marked as 'unlabeled'.

        Returns:
            A list of database rows corresponding to unlabeled structures.
        """
        with ase.db.connect(self.db_path) as db:
            return list(db.select(state="unlabeled"))

    def update_row(
        self, row_id: int, data: dict, key_value_pairs: dict
    ) -> None:
        """
        Updates a specific row in the database with new data and key-value pairs.

        Args:
            row_id: The ID of the row to update.
            data: A dictionary of results (e.g., energy, forces) to store.
            key_value_pairs: A dictionary of metadata to update.
        """
        with ase.db.connect(self.db_path) as db:
            db.update(row_id, data=data, **key_value_pairs)

    def get_labeled_atoms(self) -> list[tuple[Atoms, dict]]:
        """
        Retrieves all labeled atoms and their associated data.

        Returns:
            A list of tuples, where each tuple contains an Atoms object
            and a dictionary with its key-value pairs and DFT data.
        """
        labeled_atoms = []
        with ase.db.connect(self.db_path) as db:
            for row in db.select(state="labeled"):
                atoms = row.toatoms()
                # Combine key_value_pairs and the data column for the engine
                all_data = dict(row.key_value_pairs)
                if row.data:
                    all_data["data"] = dict(row.data)
                labeled_atoms.append((atoms, all_data))
        return labeled_atoms
