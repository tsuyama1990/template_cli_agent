"""Manages all interactions with the ASE database."""

from pathlib import Path
from typing import List
from ase.db import connect
from ase import Atoms
import logging

# Configure logger for this module
logger = logging.getLogger(__name__)

class DatabaseManager:
    """
    Handles all database operations, providing a clean interface to the ASE DB.

    This class abstracts away the specific details of the `ase.db` library,
    making the main pipeline logic cleaner and more focused on the workflow.
    """

    def __init__(self, db_path: Path) -> None:
        """
        Initializes the DatabaseManager.

        Args:
            db_path: The path to the ASE database file.
        """
        self.db_path = db_path
        self._ensure_parent_directory_exists()

    def _ensure_parent_directory_exists(self) -> None:
        """Ensures that the parent directory of the database file exists."""
        try:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            logger.exception("Error creating directory for database:")
            msg = "Error creating directory for database"
            raise OSError(msg) from e

    def write_structures(self, structures: List[Atoms]) -> None:
        """
        Writes a list of ASE Atoms objects to the database.

        Each structure is written with its associated energy, forces, and stress,
        which are expected to be present in the `atoms.info` or `atoms.arrays`
        attributes.

        Args:
            structures: A list of ASE Atoms objects to write.
        """
        if not structures:
            logger.warning("No structures provided to write to the database.")
            return

        logger.info(f"Connecting to database at: {self.db_path}")
        try:
            with connect(self.db_path) as db: # type: ignore
                for atoms in structures:
                    # ASE db automatically extracts energy, forces, stress if they
                    # are attached to the Atoms object's calculator or info dict.
                    db.write(atoms)
            logger.info(f"Successfully wrote {len(structures)} structures to the database.")
        except Exception as e:
            logger.exception("Failed to write structures to database:")
            msg = "Failed to write structures to database"
            raise RuntimeError(msg) from e
