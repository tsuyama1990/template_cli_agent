from typing import List, Type

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect


class DatabaseManager:
    """Manages writing atomic structures to an ASE database."""

    def __init__(self, db_path: str) -> None:
        self.db_path = db_path
        self.db = None

    def __enter__(self) -> "DatabaseManager":
        self.db = connect(self.db_path)  # type: ignore[no-untyped-call]
        return self

    def __exit__(
        self,
        exc_type: Type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: object | None,
    ) -> None:
        pass  # The connection is implicitly closed by ase.db

    def write_structures(self, structures: List[Atoms]) -> None:
        """Writes a list of Atoms objects to the database."""
        if self.db is None:
            msg = "Database is not connected. Use 'with' statement."
            raise ConnectionError(msg)
        for atoms in structures:
            # We add dummy energy and forces as a real calculator is not used in Cycle 1
            # Attach a SinglePointCalculator to store the results
            dummy_forces = [[0.0, 0.0, 0.0]] * len(atoms)
            calc = SinglePointCalculator(atoms, energy=0.0, forces=dummy_forces)
            atoms.calc = calc
            self.db.write(atoms)
