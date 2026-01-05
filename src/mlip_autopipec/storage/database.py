"""ASE database wrapper for storing atomic structures."""

from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

from ase.atoms import Atoms
from ase.db import connect


@contextmanager
def ase_db_wrapper(db_path: Path) -> Generator[None, None, None]:
    """
    A context manager for the ASE database.

    Args:
        db_path: Path to the database file.
    """
    db = connect(db_path)
    try:
        yield db
    finally:
        pass


def write_structures(db_path: Path, structures: list[Atoms]) -> None:
    """
    Writes a list of ASE Atoms objects to the database.

    Args:
        db_path: Path to the database file.
        structures: A list of ASE Atoms objects.
    """
    with ase_db_wrapper(db_path) as db:
        for atoms in structures:
            db.write(atoms)
