from typing import Any

import ase.db
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator


class AseDBWrapper:
    """A wrapper class for the ASE database to handle data persistence."""

    def __init__(self, db_path: str):
        """
        Initializes the AseDBWrapper.
        Args:
            db_path: The path to the ASE database file.
        """
        if not db_path:
            raise ValueError("Database path cannot be empty.")
        self.db_path = db_path

    def _connect(self) -> ase.db.core.Database:
        """Returns a connection to the database."""
        return ase.db.connect(self.db_path)

    def add_atoms(self, atoms_list: list[Atoms], **kwargs):
        """
        Adds a list of Atoms objects to the database.
        Args:
            atoms_list: A list of ASE Atoms objects.
            **kwargs: Key-value pairs to add to each row. 'labelled' will be overwritten.
        """
        kvp = kwargs.copy()
        kvp['labelled'] = False

        with self._connect() as db:
            for atoms in atoms_list:
                db.write(atoms, key_value_pairs=kvp)

    def get_row(self, row_id: int) -> ase.db.row.AtomsRow | None:
        """
        Retrieves a single AtomsRow object by its ID.
        Args:
            row_id: The ID of the row to retrieve.
        Returns:
            The AtomsRow object, or None if not found.
        """
        with self._connect() as db:
            try:
                return db.get(id=row_id)
            except KeyError:
                return None

    def get_rows_to_label(self) -> list[ase.db.row.AtomsRow]:
        """
        Retrieves rows from the database that have not been labeled yet.
        Returns:
            A list of AtomsRow objects where `labelled=False`.
        """
        with self._connect() as db:
            # Note: The ase.db.select method does not support parameterized queries.
            # This is a potential security risk if user input is ever used here.
            return list(db.select('labelled=False'))

    def update_row_with_dft_results(self, row_id: int, dft_results: dict[str, Any]):
        """
        Updates a row with DFT results and marks it as labeled.
        Args:
            row_id: The ID of the row to update.
            dft_results: A dictionary containing DFT results, must have 'energy',
                         'forces', and 'stress' keys.
        """
        with self._connect() as db:
            row = db.get(id=row_id)
            atoms = row.toatoms()

            calc = SinglePointCalculator(
                atoms,
                energy=dft_results.get('energy'),
                forces=dft_results.get('forces'),
                stress=dft_results.get('stress'),
            )
            atoms.calc = calc

            new_kvp = row.key_value_pairs
            new_kvp['labelled'] = True

            db.update(row_id, atoms=atoms, **new_kvp)


    def get_all_labeled_rows(self) -> list[ase.db.row.AtomsRow]:
        """
        Retrieves all rows that have been successfully labeled.
        Returns:
            A list of AtomsRow objects where `labelled=True`.
        """
        with self._connect() as db:
            # Note: The ase.db.select method does not support parameterized queries.
            # This is a potential security risk if user input is ever used here.
            return list(db.select('labelled=True'))
