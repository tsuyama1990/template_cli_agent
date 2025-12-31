"""Defines the abstract ports for the application's core logic."""

from typing import Protocol, List
from ase import Atoms

from .models import DFTResult


class DatabasePort(Protocol):
    """
    An abstract port for database operations, defining the interface for data persistence.
    """

    def add_atoms(self, atoms_list: List[Atoms]):
        ...

    def get_atoms_by_state(self, state: str) -> List[Atoms]:
        ...

    def update_with_dft_results(self, atoms_id: int, results: DFTResult):
        ...


class ProcessRunnerPort(Protocol):
    """
    An abstract port for running external processes.
    """

    def run(self, command: List[str], input_str: str) -> str:
        ...
