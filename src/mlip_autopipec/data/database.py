import ase.db
from ase.atoms import Atoms
from pathlib import Path
from typing import Optional, Union, Tuple, List

from .models import DFTResult

class AseDB:
    """
    A wrapper class for ASE's database functionality to provide a clear,
    type-hinted interface for writing and reading calculation results.
    """

    def __init__(self, db_path: Union[str, Path]):
        """
        Initialises the database connection.

        Args:
            db_path: The path to the SQLite database file.
        """
        self.db_path = Path(db_path)
        self._connect()

    def _connect(self):
        """Connects to the ASE database, creating the parent directory if needed."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = ase.db.connect(self.db_path)

    def write(self, atoms: Atoms, result: DFTResult) -> int:
        """
        Writes an Atoms object and its corresponding DFTResult to the database.

        The DFTResult model is converted to a dictionary. Any fields with a
        value of None are excluded to prevent errors with the ASE DB backend.

        Args:
            atoms: The ase.Atoms object representing the structure.
            result: The Pydantic DFTResult object with the calculation data.

        Returns:
            The integer ID of the newly created row in the database.
        """
        # Use model_dump to convert the Pydantic model to a dict.
        # exclude_none=True is critical as the ASE DB cannot handle None values.
        key_value_pairs = result.model_dump(exclude_none=True)
        db_id = self.connection.write(atoms, key_value_pairs=key_value_pairs)
        return db_id

    def get(self, db_id: int) -> Optional[Tuple[Atoms, DFTResult]]:
        """
        Retrieves an Atoms object and its DFTResult by its database ID.

        Args:
            db_id: The integer ID of the row to retrieve.

        Returns:
            A tuple containing the Atoms object and the reconstructed DFTResult
            object, or None if the ID is not found.
        """
        try:
            row = self.connection.get(id=db_id)
        except (KeyError, IndexError):
            return None

        atoms = row.toatoms()

        # Reconstruct the DFTResult object from the row's key-value data.
        # Using .get() provides default None values for missing keys, which
        # is safe for the optional fields in the Pydantic model.
        dft_result_data = {
            'total_energy_ev': row.get('total_energy_ev'),
            'forces': row.get('forces'),
            'stress': row.get('stress'),
            'was_successful': row.get('was_successful'),
            'error_message': row.get('error_message'),
        }

        # Filter out keys that were not present in the database row to avoid
        # passing None for fields that don't have a default.
        # This makes the reconstruction robust.
        reconstruct_data = {k: v for k, v in dft_result_data.items() if k in row}

        result = DFTResult(**reconstruct_data)

        return atoms, result

    def get_all_atoms(self) -> List[Atoms]:
        """
        Retrieves all Atoms objects from the database.

        Returns:
            A list of all ase.Atoms objects stored in the database.
        """
        return [row.toatoms() for row in self.connection.select()]
