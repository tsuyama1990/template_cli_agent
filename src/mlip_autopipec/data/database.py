from typing import List, Tuple
import ase.db
from ase import Atoms

from mlip_autopipec.data.models import DFTResult


class AseDB:
    """
    A wrapper class for the ASE database to provide a high-level API for
    database operations. This class uses a connect-on-demand pattern
    to ensure transaction integrity.
    """

    def __init__(self, db_path: str):
        """
        Initializes the AseDB wrapper with the path to the database.

        Args:
            db_path: The path to the SQLite database file.
        """
        self.db_path = db_path

    def add_atoms(self, atoms: Atoms, state: str = "initial") -> int:
        """
        Adds a new ASE Atoms object to the database.

        Args:
            atoms: The ASE Atoms object to add.
            state: The initial state of the structure.

        Returns:
            The unique integer ID assigned to the new row.
        """
        with ase.db.connect(self.db_path) as con:
            return con.write(atoms, key_value_pairs={"state": state})

    def get_atoms_to_label(self, limit: int = 10) -> List[Atoms]:
        """
        Fetches a list of structures from the database that need to be
        processed by the Labelling Engine.

        Args:
            limit: The maximum number of structures to fetch.

        Returns:
            A list of ASE Atoms objects.
        """
        with ase.db.connect(self.db_path) as con:
            rows = list(con.select(state="initial", limit=limit))
        return [row.toatoms() for row in rows]

    def write_dft_result(self, atoms_id: int, result: DFTResult):
        """
        Writes the complete outcome of a DFT calculation to the database.

        Args:
            atoms_id: The ID of the atoms row to update.
            result: The DFTResult object to write.
        """
        dft_result_json = result.model_dump_json()
        with ase.db.connect(self.db_path) as con:
            con.update(atoms_id, state="labelled", dft_result=dft_result_json)

    def update_state(self, atoms_id: int, new_state: str):
        """
        Provides a generic mechanism to update the state of a structure.

        Args:
            atoms_id: The ID of the atoms row to update.
            new_state: The new state to set.
        """
        with ase.db.connect(self.db_path) as con:
            con.update(atoms_id, state=new_state)

    def get_training_data(self) -> Tuple[List[Atoms], List[DFTResult]]:
        """
        Retrieves a complete dataset of all successfully labelled atoms
        and their corresponding results.

        Returns:
            A tuple containing a list of ASE Atoms objects and a list of
            DFTResult objects.
        """
        with ase.db.connect(self.db_path) as con:
            rows = list(con.select(state="labelled"))
        atoms_list = [row.toatoms() for row in rows]
        results_list = [
            DFTResult.model_validate_json(row.key_value_pairs["dft_result"])
            for row in rows
        ]
        return atoms_list, results_list
