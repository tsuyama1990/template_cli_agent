from typing import Any

from ase.atoms import Atoms
from ase.db import connect


class AseDB:
    """
    A wrapper class for the ASE database to handle storing and retrieving atomic data.
    """

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.connection = connect(db_path)

    def write(self, atoms: Atoms, dft_result: dict[str, Any]) -> int:
        """
        Writes an Atoms object and its corresponding DFT results to the database.

        Args:
            atoms: The ASE Atoms object.
            dft_result: A dictionary containing the DFT calculation results.

        Returns:
            The ID of the new row in the database.
        """
        # ASE db cannot handle None values for key_value_pairs
        filtered_dft_result = {k: v for k, v in dft_result.items() if v is not None}
        return self.connection.write(atoms, key_value_pairs=filtered_dft_result)

    def get(self, selection) -> tuple[Atoms, dict[str, Any]]:
        """
        Retrieves an entry from the database.

        Args:
            selection: The selection criteria for the database query.

        Returns:
            A tuple containing the Atoms object and its key-value pairs.
        """
        row = self.connection.get(selection)
        if row:
            return row.toatoms(), row.key_value_pairs
        return None, None
