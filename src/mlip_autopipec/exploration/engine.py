"""Placeholder for the MD Exploration Engine."""
from typing import List

from ase.atoms import Atoms


class MDExplorer:
    """A placeholder class for the MD exploration engine."""

    def run_md(self, structures: List[Atoms]) -> List[Atoms]:
        """
        A pass-through method that returns the input structures without modification.

        Args:
            structures: A list of ase.Atoms objects.

        Returns:
            The same list of ase.Atoms objects.
        """
        print("Executing placeholder MD exploration...")
        return structures
