from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


class MDExplorer:
    """Placeholder for the MD Explorer component."""

    def run_md(self, structures: list[Atoms]) -> list[Atoms]:
        """A pass-through method that returns the input structures."""
        print("Skipping MD exploration in Cycle 1.")
        return structures
