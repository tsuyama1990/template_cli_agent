from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.db.core import Database
from ase.stress import voigt_6_to_full_3x3_stress

from mlip_autopipec.domain_models import DFTResult


class AseDB:
    def __init__(self, db_path: Path):
        self.db_path = db_path

    @contextmanager
    def get_connection(self) -> Generator[Database, None, None]:
        """Provides a managed connection to the ASE database."""
        db = connect(self.db_path)
        try:
            yield db
        finally:
            pass  # The connection is automatically handled by ase.db.connect

    def write(
        self,
        atoms: Atoms,
        dft_result: DFTResult,
        was_successful: bool,
        key_value_pairs: dict[str, Any] | None = None,
    ) -> int:
        """Writes an Atoms object and its calculation results to the database."""
        if dft_result.stress is not None:
            # ASE calculators expect stress in Voigt format
            if len(dft_result.stress) == 3:  # Full 3x3 matrix
                stress_voigt = dft_result.stress
            else:
                stress_voigt = voigt_6_to_full_3x3_stress(dft_result.stress)
        else:
            stress_voigt = None

        calculator = SinglePointCalculator(
            atoms,
            energy=dft_result.energy,
            forces=dft_result.forces,
            stress=stress_voigt,
        )
        atoms.calc = calculator

        with self.get_connection() as db:
            uid = db.write(
                atoms,
                key_value_pairs={
                    "was_successful": was_successful,
                    **(key_value_pairs or {}),
                },
            )
        return uid

    def get(self, uid: int) -> Atoms | None:
        """Retrieves an Atoms object from the database by its unique ID."""
        with self.get_connection() as db:
            try:
                row = db.get(id=uid)
                return row.toatoms()
            except KeyError:
                return None
