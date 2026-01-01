import ase.db
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import full_3x3_to_voigt_6_stress

from mlip_autopipec.data.models import DFTResults


class AseDBWrapper:
    """
    A wrapper class for the ASE database to handle data persistence.

    This wrapper ensures that all interactions with the ASE database are
    handled consistently and safely.
    """

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
        """Returns a new connection to the database."""
        return ase.db.connect(self.db_path)

    def add_structures(self, atoms_list: list[Atoms], **kwargs):
        """
        Adds a list of Atoms objects (structures) to the database.

        Args:
            atoms_list: A list of ASE Atoms objects.
            **kwargs: Key-value pairs to add to each row. The 'labelled' key
                      will be overwritten.
        """
        kvp = kwargs.copy()
        kvp['labelled'] = False

        with self._connect() as db:
            for atoms in atoms_list:
                db.write(atoms, key_value_pairs=kvp)

    def get_row(self, row_id: int) -> ase.db.row.AtomsRow | None:
        """
        Retrieves a single AtomsRow object by its ID.
        """
        with self._connect() as db:
            try:
                return db.get(id=row_id)
            except KeyError:
                return None

    def get_rows_to_label(self) -> list[ase.db.row.AtomsRow]:
        """
        Retrieves rows from the database that have not been labeled yet.
        """
        with self._connect() as db:
            return list(db.select(labelled=False))

    def update_row_with_dft_results(self, row_id: int, dft_results: DFTResults):
        """
        Updates a row with DFT results and marks it as labeled.
        """
        with self._connect() as db:
            row = db.get(id=row_id)
            atoms = row.toatoms()

            # The SinglePointCalculator requires stress in 6-element Voigt form.
            voigt_stress = full_3x3_to_voigt_6_stress(dft_results.stress)

            calc = SinglePointCalculator(
                atoms,
                energy=dft_results.energy,
                forces=dft_results.forces,
                stress=voigt_stress,
            )
            atoms.calc = calc

            new_kvp = row.key_value_pairs
            new_kvp['labelled'] = True

            db.update(row_id, atoms=atoms, **new_kvp)

    def get_all_labeled_rows(self) -> list[ase.db.row.AtomsRow]:
        """
        Retrieves all rows that have been successfully labeled.
        """
        with self._connect() as db:
            return list(db.select(labelled=True))
