# src/mlip_autopipec/interfaces.py
from typing import Any, Protocol

from ase import Atoms


class IStructureGenerator(Protocol):
    """Interface for structure generators."""

    def generate(self) -> list[Atoms]:
        """Generate a list of atomic structures."""
        ...


class ISampler(Protocol):
    """Interface for samplers."""

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """Sample structures from a trajectory."""
        ...


class IDatabaseManager(Protocol):
    """Interface for the database manager."""

    def connect(self, path: str) -> Any:
        """Connect to the database."""
        ...

    def write_structures(self, structures: list[Atoms]) -> None:
        """Write structures to the database."""
        ...


class IPipelineOrchestrator(Protocol):
    """Interface for the pipeline orchestrator."""

    def run(self) -> None:
        """Run the pipeline."""
        ...
